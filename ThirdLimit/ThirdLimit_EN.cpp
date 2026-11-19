// https://www.boost.org/
// used for high-precision calculations
// I have boost_1_88_0-msvc-14.3-64.exe

#include <vector>
#include <iostream>
#include <fstream>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/statistics/linear_regression.hpp>
#include <boost/format.hpp>
#include <Windows.h>
#include <tuple>

using namespace boost::multiprecision;
using namespace std;

cpp_dec_float_50 pi = 3.14159265358979323846264338327950288419716939937510;

typedef tuple <cpp_dec_float_100, cpp_dec_float_100, cpp_dec_float_100> ids;

int main() {

// We're looking at 7 points. We could probably use 5, or even 15.
    std::vector<cpp_dec_float_100> x(7);       // number of angles; square, octagon, sixteen, thirty-two, 64, 128 ...
    std::vector<cpp_dec_float_100> alfa(7);    // alpha angle
    std::vector<cpp_dec_float_100> alfaLog(7); // logarithm of the alpha angle. Comes in handy.
    std::vector<cpp_dec_float_100> beta(7);    // beta angle
    std::vector<cpp_dec_float_100> exp(7);     // power of two

    std::ofstream outFile;
    outFile.open("d:\\temp\\ThirdLimitOut.txt"); // d:\temp must exist

    cpp_dec_float_100 n;
    cpp_dec_float_100 m;

    outFile << "vertices A, B, C of triangle ABC tend to point C.\n";
    outFile << "thus, all three points A, B, C become point ABC or simply C\n";
    outFile << "maybe.\n";

    outFile << "(1) Limit lim(360 / n) = 0 at n->8 achieved, but when I did not find out. triangles AiOC\n";
    outFile << "(2) Limit lim(180 - 360 / n) = 180 at n->8 is reached at the 332nd point.. triangles AiBiC\n";
    outFile << std::endl;

    m = 331;
    n = pow(2, m);
    outFile << "power of two: m = " << m << '\n';
    outFile << std::setprecision(std::numeric_limits<cpp_dec_float_100>::digits10) << "n = pow(2, m); n = " << n << " angles\n";
    outFile << "beta = 180. - 360. / n;\n";
    outFile << "beta = " << 180. - 360. / n << " - this is still a triangle АВС,\n";
    outFile << "                                                                                                         very small and very obtuse\n";
    outFile << std::endl;

    m = 332;
    n = pow(2, m);
    outFile << "power of two: m = " << m << '\n';
    outFile << std::setprecision(std::numeric_limits<cpp_dec_float_100>::digits10) << "n = pow(2, m); n = " << n << '\n';
    outFile << "beta = 180. - 360. / n;\n";
    outFile << "beta = " << 180. - 360. / n << " - reached the limit. Triangle ABC apparently became point C.\n";
    outFile << "perhaps this is just a computational quirk, but here's the quirk. After all, 2 to the power of 332 is a polygon, which is a lot of angles.)\n";

    outFile << std::endl;

    n = 4; // number of corners starting from the square
    int ind = 0;

    for (int i = 2; i < 9; i++)
    {
        x[ind] = n; // number of angles
        alfa[ind] = 360. / n;
        alfaLog[ind] = std::log10((double)alfa[ind]);
        beta[ind] = 180. - alfa[ind];
        exp[ind] = i;

        ind++;
        n *= 2; // = pow(2, i); // you can (should) try with 3 and then 6, 12, 24, etc.
    }

    outFile << "exp  n      alfa        beta  alfaLog - для _least_squares_with_R_squared\n";

    for (int i = 0; i < 7; i++)
    {
        outFile << exp[i] << "  " << std::setw(3) << x[i] << "  " << std::setw(8) << alfa[i] << "  "
            << std::setw(10) << beta[i] << "  " << std::setw(10) << alfaLog[i] << std::endl;
    }
    outFile << std::endl;

    ids c0c1R2; // <-everything at once for f(x) = c0 + c1 * x и R2
    outFile << "attempt to study sequences (1) и (2)\n" << std::endl;

    outFile << "f(x) = c0 + c1 * x" << std::endl;
    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(exp, alfaLog);
    outFile << "_least_squares_with_R2(exp, alfaLog);" << std::endl;
    outFile << "c0 = " << get<0>(c0c1R2) << "\n"; outFile << "c1 = " << get<1>(c0c1R2) << "\n"; outFile << "R2 = " << get<2>(c0c1R2) << " - good exponent\n";
    outFile << std::endl;

    outFile << "We are interested in the side of the polygon, so we will take beta i+1 to alpha i." << std::endl;
    outFile << std::endl;

    n = 4; // number of corners starting from the square
    for (int i = 0; i < 7; i++)
    {
        n *= 2; // beta i+1 to alpha i - this is here, first doubled, then calculated
        beta[i] = 180. - 360. / n;
    }

    std::vector<cpp_dec_float_100> xp(7);
    for (int i = 0; i < 7; i++)
    {
        xp[i] = 4 / x[i]; // I found this scale by accident.
    }

    outFile << "xp[i] = 4 / x[i];" << std::endl;
    outFile << "i          xp        alfa  beta" << std::endl;
    for (int i = 0; i < 7; i++)
    {
        outFile << i << "  " << std::setw(10) << xp[i] << "  " << std::setw(10) << alfa[i] << "  " << beta[i] << std::endl;
    }
    outFile << std::endl;

    outFile << "I accidentally found the scale xp.\n";
    outFile << "_least_squares_with_R2(xp, alfa);" << std::endl;
    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(xp, alfa);
    outFile << "c0 = " << get<0>(c0c1R2) << "    ";       outFile << "c1 = " << get<1>(c0c1R2) << "    "; outFile << "R2 = " << get<2>(c0c1R2) << " (1)\n";
    outFile << "_least_squares_with_R2(xp, beta);" << std::endl;
    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(xp, beta);
    outFile << "c0 = " << get<0>(c0c1R2) << "          "; outFile << "c1 = " << get<1>(c0c1R2) << "   ";  outFile << "R2 = " << get<2>(c0c1R2) << " (2)\n";

    outFile << std::endl;

    outFile << "Note that limits (1) and (2) behave differently.." << std::endl;
    outFile << std::endl;

    ind = 0;
    for (int i = 2; ind < 7; ind++)
    {
        exp[ind] = i;
        n = pow(2, i);
        alfa[ind] = 360. / n;
        xp[ind] = 4 / n;
        n = pow(2, i + 1); // for alpha and beta i+1
        beta[ind] = 180. - 360. / n;
        i += 10;
    }

    outFile << "n = pow(2, i); i = "; for (auto& c : exp) outFile << c << "  "; outFile << "with a step i += 10\n";

    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(xp, alfa);
    outFile << "c0 = " << std::setw(8) << get<0>(c0c1R2) << "    "; outFile << "c1 = " << std::setw(3) << get<1>(c0c1R2) << "    "; outFile << "R2 = " << get<2>(c0c1R2) << "  (1)\n";
    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(xp, beta);
    outFile << "c0 = " << std::setw(8) << get<0>(c0c1R2) << "    "; outFile << "c1 = " << std::setw(3) << get<1>(c0c1R2) << "    "; outFile << "R2 = " << get<2>(c0c1R2) << "  (2)\n";
    outFile << std::endl;

    ind = 0;
    for (int i = 2; ind < 7; ind++)
    {
        exp[ind] = i;
        n = pow(2, i);
        alfa[ind] = 360. / n;
        xp[ind] = 4 / n;
        n = pow(2, i + 1);
        beta[ind] = 180. - 360. / n;
        i += 25;
    }
    outFile << "n = pow(2, i); i = "; for (auto& c : exp) outFile << c << "  "; outFile << "i += 25\n";

    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(xp, alfa);
    outFile << "c0 = " << std::setw(8) << get<0>(c0c1R2) << "    "; outFile << "c1 = " << std::setw(3) << get<1>(c0c1R2) << "    "; outFile << "R2 = " << get<2>(c0c1R2) << "  (1)\n";
    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(xp, beta);
    outFile << "c0 = " << std::setw(8) << get<0>(c0c1R2) << "    "; outFile << "c1 = " << std::setw(3) << get<1>(c0c1R2) << "    "; outFile << "R2 = " << get<2>(c0c1R2) << "  (2)\n";
    outFile << std::endl;

    for (int i = 2, ind = 0; ind < 7; ind++)
    {
        exp[ind] = i;
        n = pow(2, i);
        alfa[ind] = 360. / n;
        xp[ind] = 4 / n;
        n = pow(2, i + 1);
        beta[ind] = 180. - 360. / n;
        i += 48;
    }
    outFile << "n = pow(2, i); i = "; for (auto& c : exp) outFile << c << "  "; outFile << "i += 48\n";

    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(xp, alfa);
    outFile << "c0 = " << std::setw(8) << get<0>(c0c1R2) << "    "; outFile << "c1 = " << std::setw(3) << get<1>(c0c1R2) << "    "; outFile << "R2 = " << get<2>(c0c1R2) << "  (1)\n";
    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(xp, beta);
    outFile << "c0 = " << std::setw(8) << get<0>(c0c1R2) << "    "; outFile << "c1 = " << std::setw(3) << get<1>(c0c1R2) << "    "; outFile << "R2 = " << get<2>(c0c1R2) << "  (2)\n";
    outFile << std::endl;

    outFile << "--------------------------------------------------------------------------------------------------------------------------------\n\n";
    outFile << "Side length" << std::endl;
    outFile << std::endl;
    outFile << "for verification" << std::endl;

    outFile << "2 * boost::multiprecision::sin(pi / 6) inscribed hexagon. must be = 1" << std::endl;
    outFile << 2 * boost::multiprecision::sin(pi / 6) << " - аlmost one, but not one. It's probably because of the sin." << std::endl << std::endl;

    boost::multiprecision::cpp_bin_float_100 x2 = 2;
    //outFile << boost::multiprecision::sqrt(2) << " - to check 2, just 2" << std::endl; - You can't do that. It calculates incorrectly. The result is 1.
    outFile << "boost::multiprecision::sqrt(2) = " << std::endl;
    outFile << "= " << boost::multiprecision::sqrt(x2) << std::endl;
    outFile << "  1,4142135623730950488 - from the Internet\n";

    outFile << std::endl;

    n = 4; // number of corners starting from the square
    ind = 0;
    std::vector<cpp_dec_float_100> a(7);       // side length
    std::vector<cpp_dec_float_100> aLog(7);       // log side length

    outFile << "   n  2 * sin(pi / n)" << std::endl;
    for (int i = 2; i < 9; i++)
    {
        a[ind] = 2 * boost::multiprecision::sin(pi / n);
        aLog[ind] = log10(a[ind]);
        exp[ind] = i;

        outFile << setw(4) << n << "  " << a[ind] << std::endl;

        ind++;
        n *= 2;
    }

    outFile << std::endl;

    outFile << "_least_squares_with_R_squared(xp, a);" << std::endl;
    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(xp, a);
    outFile << "c0 = " << get<0>(c0c1R2) << "\n"; outFile << "c1 = " << get<1>(c0c1R2) << "\n"; outFile << "R2 = " << get<2>(c0c1R2) << "\n";

    outFile << std::endl;

    outFile << "_least_squares_with_R_squared(exp, aLog);" << std::endl;
    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(exp, aLog);
    outFile << "c0 = " << get<0>(c0c1R2) << "\n"; outFile << "c1 = " << get<1>(c0c1R2) << "\n"; outFile << "R2 = " << get<2>(c0c1R2) << "\n";

    outFile << std::endl;

    cpp_dec_float_100 rd;
    cpp_dec_float_100 c;

    outFile << "for verification. hexagon inscribed in a circle" << std::endl;

    outFile << "C = 2 * sin(pi / 6) =         " << 2 * sin(pi / 6) << " side of a 6-gon inscribed in a circle\n";
    rd = 360 / 6 * pi / 180;
    c = 2 * sin(rd / 2);
    outFile << "C = 2 * sin(360 / 6) =        " << c << " side of an isosceles triangle with an angle of 60 degrees\n";
    rd = 360 / 12 * pi / 180;
    c = 2 * sin(rd / 2);
    rd = (180 - 360 / 12) * pi / 180;
    outFile << "C = 2 * sin(180 - 360 / 12) = " << 2 * c * sin(rd / 2) << " side of an isosceles triangle with an angle on the circle\n";

    outFile << std::endl;
    outFile << "number of angles       side length" << std::endl;

    int step = 5;
    for (int i = 0; i < 10; i++)
    {
        int e = 2 + i * step;
        rd = 360 / pow(2, e) * pi / 180; // alpha angle in radians
        c = 2 * sin(rd / 2);             // side C from side alpha
        outFile << setw(15) << pow(2, e) << "  " << c << "  alfa\n";
        e++;
        rd = 360 / pow(2, e) * pi / 180; // next alpha angle in radians
        c = 2 * sin(rd / 2);             // and side C for calculating C from side beta
        rd = (180 - 360 / pow(2, e)) * pi / 180; // angle beta in radians
        outFile << setw(15) << pow(2, e) << "  " << 2 * c * sin(rd / 2) << "  beta\n";
    }
    outFile << std::endl;

    outFile << "-------------------------" << std::endl;
    outFile << std::endl;

    outFile << "point size" << std::endl;
    outFile << std::endl;

    outFile << "C = 2 * sin(pi / pow(2, 332)) =  " << 2 * sin(pi / pow(2, 332)) << " side of a polygon inscribed in a circle\n";
    rd = 360 / pow(2, 332) * pi / 180;
    c = 2 * sin(rd / 2);
    outFile << "C = 2 * sin(360 / pow(2, 332)) = " << c << " side of an isosceles triangle with angle alpha\n";

    rd = 360 / pow(2, 333) * pi / 180;
    c = 2 * sin(rd / 2);
    rd = (180 - 360 / pow(2, 333)) * pi / 180;
    outFile << "C = 2 * c * sin(180 - beta) =    " << 2 * c * sin(rd / 2) << " side of an isosceles triangle with angle beta\n";
    outFile << std::endl;

    std::cout << boost::format("\nC = 2 * c * sin(180 - beta) = %s\n") % (2 * c * sin(rd / 2));

    outFile << "I'm here..";
    
    outFile.close();

    return 2;

// translated by Google
}