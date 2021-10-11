#include "BezierEasing.h"

// These values are established by empiricism with tests (tradeoff: performance VS precision)
static const int NEWTON_ITERATIONS = 4;
static const double NEWTON_MIN_SLOPE = 0.001;
static const double SUBDIVISION_PRECISION = 0.0000001;
static const int SUBDIVISION_MAX_ITERATIONS = 10;

static const int kSplineTableSize = 11;
static const double kSampleStepSize = 1.0 / ((double)kSplineTableSize - 1.0);

static double A(double aA1, double aA2) { return 1.0 - 3.0 * aA2 + 3.0 * aA1; }
static double B(double aA1, double aA2) { return 3.0 * aA2 - 6.0 * aA1; }
static double C(double aA1) { return 3.0 * aA1; }

// Returns x(t) given t, x1, and x2, or y(t) given t, y1, and y2.
static double calcBezier(double aT, double aA1, double aA2) { return ((A(aA1, aA2) * aT + B(aA1, aA2)) * aT + C(aA1)) * aT; }

// Returns dx/dt given t, x1, and x2, or dy/dt given t, y1, and y2.
static double getSlope(double aT, double aA1, double aA2) { return 3.0 * A(aA1, aA2) * aT * aT + 2.0 * B(aA1, aA2) * aT + C(aA1); }

static double binarySubdivide(double aX, double aA, double aB, double mX1, double mX2) {
    double currentX, currentT;
    int i = 0;

    do {
        currentT = aA + (aB - aA) / 2.0;
        currentX = calcBezier(currentT, mX1, mX2) - aX;

        if (currentX > 0.0)
            aB = currentT;
        else
            aA = currentT;
    } while (abs(currentX) > SUBDIVISION_PRECISION && ++i < SUBDIVISION_MAX_ITERATIONS);

    return currentT;
}

static double newtonRaphsonIterate(double aX, double aGuessT, double mX1, double mX2) {
    for (int i = 0; i < NEWTON_ITERATIONS; ++i) {
        double currentSlope = getSlope(aGuessT, mX1, mX2);

        if (currentSlope == 0.0)
            return aGuessT;

        double currentX = calcBezier(aGuessT, mX1, mX2) - aX;
        aGuessT -= currentX / currentSlope;
    }

    return aGuessT;
}

static double LinearEasing(double x) {
    return x;
}

std::function<double(double)> BezierEasing(double mX1, double mY1, double mX2, double mY2) {
    if (!(0 <= mX1 && mX1 <= 1 && 0 <= mX2 && mX2 <= 1))
        throw "bezier x values must be in[0, 1] range";

    if (mX1 == mY1 && mX2 == mY2)
        return LinearEasing;

    // Precompute samples table
    static double sampleValues[kSplineTableSize];

    for (int i = 0; i < kSplineTableSize; ++i)
        sampleValues[i] = calcBezier(i * kSampleStepSize, mX1, mX2);

    auto getTForX = [=](double aX) -> double {
        double intervalStart = 0.0;
        int currentSample = 1;
        int lastSample = kSplineTableSize - 1;

        for (; currentSample != lastSample && sampleValues[currentSample] <= aX; ++currentSample)
            intervalStart += kSampleStepSize;

        --currentSample;

        // Interpolate to provide an initial guess for t
        double dist = (aX - sampleValues[currentSample]) / (sampleValues[currentSample + 1] - sampleValues[currentSample]);
        double guessForT = intervalStart + dist * kSampleStepSize;
        double initialSlope = getSlope(guessForT, mX1, mX2);

        if (initialSlope >= NEWTON_MIN_SLOPE)
            return newtonRaphsonIterate(aX, guessForT, mX1, mX2);
        else if (initialSlope == 0.0)
            return guessForT;
        else
            return binarySubdivide(aX, intervalStart, intervalStart + kSampleStepSize, mX1, mX2);
    };

    return [=](double x) -> double {
        // Because JavaScript number are imprecise, we should guarantee the extremes are right.
        if (x == 0 || x == 1)
            return x;

        return calcBezier(getTForX(x), mY1, mY2);
    };
};