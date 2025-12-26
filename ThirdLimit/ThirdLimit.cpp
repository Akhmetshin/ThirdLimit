// https://www.boost.org/
// использую для вычислений с повышенной точностью
// у меня boost_1_88_0-msvc-14.3-64.exe

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

    // смотрим 7 точек. наверное, можно 5, а можно и 15.
    std::vector<cpp_dec_float_100> x(7);       // количество углов; квадрат, восьми[угольник], шестнадцати[], тридцатидвух[], 64, 128 ...
    std::vector<cpp_dec_float_100> alfa(7);    // угол альфа
    std::vector<cpp_dec_float_100> alfaLog(7); // логарифм угла альфа. пригодится.
    std::vector<cpp_dec_float_100> beta(7);    // угол бетта
    std::vector<cpp_dec_float_100> exp(7);     // степень двойки

    std::ofstream outFile;
    outFile.open("d:\\temp\\ThirdLimitOut.txt"); // d:\temp должна существовать

    cpp_dec_float_100 n;
    cpp_dec_float_100 m;

    outFile << "вершины А, В, С треугольника АВС стремятся к точке С.\n";
    outFile << "таким образом, все три точки А, В, С становятся точкой АВС или просто С\n";
    outFile << "наверное.\n";

    outFile << "(1) Предел lim(360 / n) = 0 при n->8 достигается, но когда я не выяснил. треугольники AiOC\n";
    outFile << "(2) Предел lim(180 - 360 / n) = 180 при n->8 достигается на 332-рой точке.. треугольники AiBiC\n";
    outFile << std::endl;

    m = 331;
    n = pow(2, m);
    outFile << "степень двойки: m = " << m << '\n';
    outFile << std::setprecision(std::numeric_limits<cpp_dec_float_100>::digits10) << "n = pow(2, m); n = " << n << " углов\n";
    outFile << "beta = 180. - 360. / n;\n";
    outFile << "beta = " << 180. - 360. / n << " - это ещё треугольник АВС,\n";
    outFile << "                                                                                                         очень маленький и очень тупой\n";
    outFile << std::endl;

    m = 332;
    n = pow(2, m);
    outFile << "степень двойки: m = " << m << '\n';
    outFile << std::setprecision(std::numeric_limits<cpp_dec_float_100>::digits10) << "n = pow(2, m); n = " << n << '\n';
    outFile << "beta = 180. - 360. / n;\n";
    outFile << "beta = " << 180. - 360. / n << " - достигли предела. треугольник АВС, видимо, стал точкой С\n";
    outFile << "возможно, это просто вычислительный казус, но вот такой казус. Всё-таки, 2 в степени 332 - угольник, это очень много углов.)\n";

    outFile << std::endl;

    n = 4; // количество углов начиная с квадрата
    int ind = 0;

    for (int i = 2; i < 9; i++)
    {
        x[ind] = n; // количество углов
        alfa[ind] = 360. / n;
        alfaLog[ind] = std::log10((double)alfa[ind]);
        beta[ind] = 180. - alfa[ind];
        exp[ind] = i;

        ind++;
        n *= 2; // = pow(2, i); // можно (нужно) попробовать с 3 и далее 6, 12, 24 итд.
    }

    outFile << "exp  n      alfa        beta  alfaLog - для _least_squares_with_R_squared\n";

    for (int i = 0; i < 7; i++)
    {
        outFile << exp[i] << "  " << std::setw(3) << x[i] << "  " << std::setw(8) << alfa[i] << "  "
            << std::setw(10) << beta[i] << "  " << std::setw(10) << alfaLog[i] << std::endl;
    }
    outFile << std::endl;

    ids c0c1R2; // <-сразу всё для f(x) = c0 + c1 * x и R2
    outFile << "попытка исследования последовательностей (1) и (2)\n" << std::endl;

    outFile << "f(x) = c0 + c1 * x" << std::endl;
    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(exp, alfaLog);
    outFile << "_least_squares_with_R2(exp, alfaLog);" << std::endl;
    outFile << "c0 = " << get<0>(c0c1R2) << "\n"; outFile << "c1 = " << get<1>(c0c1R2) << "\n"; outFile << "R2 = " << get<2>(c0c1R2) << " - хорошая экспонента\n";
    outFile << std::endl;

    outFile << "Нас интересует сторона многоугольника, поэтому будем брать beta i+1 к alfa i." << std::endl;
    outFile << std::endl;

    n = 4; // количество углов начиная с квадрата
    for (int i = 0; i < 7; i++)
    {
        n *= 2; // beta i+1 к alfa i - это здесь, сначала удвоили, потом рассчитали
        beta[i] = 180. - 360. / n;
    }

    std::vector<cpp_dec_float_100> xp(7);
    for (int i = 0; i < 7; i++)
    {
        xp[i] = 4 / x[i]; // случайно нашёл эту шкалу.
    }

    outFile << "xp[i] = 4 / x[i];" << std::endl;
    outFile << "i          xp        alfa  beta" << std::endl;
    for (int i = 0; i < 7; i++)
    {
        outFile << i << "  " << std::setw(10) << xp[i] << "  " << std::setw(10) << alfa[i] << "  " << beta[i] << std::endl;
    }
    outFile << std::endl;

    outFile << "случайно нашёл шкалу xp.\n";
    outFile << "_least_squares_with_R2(xp, alfa);" << std::endl;
    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(xp, alfa);
    outFile << "c0 = " << get<0>(c0c1R2) << "    ";       outFile << "c1 = " << get<1>(c0c1R2) << "    "; outFile << "R2 = " << get<2>(c0c1R2) << " (1)\n";
    outFile << "_least_squares_with_R2(xp, beta);" << std::endl;
    c0c1R2 = boost::math::statistics::simple_ordinary_least_squares_with_R_squared(xp, beta);
    outFile << "c0 = " << get<0>(c0c1R2) << "          "; outFile << "c1 = " << get<1>(c0c1R2) << "   ";  outFile << "R2 = " << get<2>(c0c1R2) << " (2)\n";

    outFile << std::endl;

    outFile << "заметим, что пределы (1) и (2) ведут себя по-разному." << std::endl;
    outFile << std::endl;

    ind = 0;
    for (int i = 2; ind < 7; ind++)
    {
        exp[ind] = i;
        n = pow(2, i);
        alfa[ind] = 360. / n;
        xp[ind] = 4 / n;
        n = pow(2, i + 1); // для alfa i  beta i+1
        beta[ind] = 180. - 360. / n;
        i += 10;
    }

    outFile << "n = pow(2, i); i = "; for (auto& c : exp) outFile << c << "  "; outFile << "с шагом i += 10\n";

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
    outFile << "Длина стороны" << std::endl;
    outFile << std::endl;
    outFile << "для проверки" << std::endl;

    outFile << "2 * boost::multiprecision::sin(pi / 6) вписанный шестиугольник. д.б. = 1" << std::endl;
    outFile << 2 * boost::multiprecision::sin(pi / 6) << " - почти единица, но не единица. наверное, это из-за синуса" << std::endl << std::endl;

    boost::multiprecision::cpp_bin_float_100 x2 = 2;
    //outFile << boost::multiprecision::sqrt(2) << " - для проверки 2, просто 2" << std::endl; - так делать нельзя. вычисляет неправильно. результат =1
    outFile << "boost::multiprecision::sqrt(2) = " << std::endl;
    outFile << "= " << boost::multiprecision::sqrt(x2) << std::endl;
    outFile << "  1,4142135623730950488 - из интернета\n";

    outFile << std::endl;

    n = 4; // количество углов начиная с квадрата
    ind = 0;
    std::vector<cpp_dec_float_100> a(7);       // длина стороны
    std::vector<cpp_dec_float_100> aLog(7);       // log длина стороны

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

    outFile << "для проверки. шестиугольник вписанный в окружность" << std::endl;

    outFile << "C = 2 * sin(pi / 6) =         " << 2 * sin(pi / 6) << " сторона 6-угольника вписанного в окружность\n";
    rd = 360 / 6 * pi / 180;
    c = 2 * sin(rd / 2);
    outFile << "C = 2 * sin(360 / 6) =        " << c << " сторона равнобедренного треугольника с углом 60 гр\n";
    rd = 360 / 12 * pi / 180;
    c = 2 * sin(rd / 2);
    rd = (180 - 360 / 12) * pi / 180;
    outFile << "C = 2 * sin(180 - 360 / 12) = " << 2 * c * sin(rd / 2) << " сторона равнобедренного треугольника с углом на окружности\n";

    outFile << std::endl;
    outFile << "колич углов       длина стороны" << std::endl;

    int step = 5;
    for (int i = 0; i < 10; i++)
    {
        int e = 2 + i * step;
        rd = 360 / pow(2, e) * pi / 180; // угол альфа в радианах
        c = 2 * sin(rd / 2);             // сторона С со стороны альфа
        outFile << setw(15) << pow(2, e) << "  " << c << "  alfa\n";
        e++;
        rd = 360 / pow(2, e) * pi / 180; // следующий угол альфа в радианах
        c = 2 * sin(rd / 2);             // и сторона С для рассчёта С со стороны бета
        rd = (180 - 360 / pow(2, e)) * pi / 180; // угол бета в радианах
        outFile << setw(15) << pow(2, e) << "  " << 2 * c * sin(rd / 2) << "  beta\n";
    }
    outFile << std::endl;

    outFile << "-------------------------" << std::endl;
    outFile << std::endl;

    outFile << "размер точки" << std::endl;
    outFile << std::endl;

    outFile << "C = 2 * sin(pi / pow(2, 332)) =            " << 2 * sin(pi / pow(2, 332)) << " сторона многоугольника вписанного в окружность\n";
    rd = 360 / pow(2, 332) * pi / 180;
    c = 2 * sin(rd / 2);
    outFile << "C = 2 * sin(360 / pow(2, 332)) =           " << c << " сторона равнобедренного треугольника с углом альфа\n";

    rd = 360 / pow(2, 333) * pi / 180;
    c = 2 * sin(rd / 2);
    rd = (180 - 360 / pow(2, 333)) * pi / 180;
    outFile << "C = 2 * c * sin(180 - beta) =              " << 2 * c * sin(rd / 2) << " сторона равнобедренного треугольника с углом бета\n";
    outFile << std::endl;

    std::cout << boost::format("\nC = 2 * c * sin(180 - 360 / pow(2, 333)) = %s\n") % (2 * c * sin(rd / 2));

    outFile << "я здесь..";
    
    return 2;

}