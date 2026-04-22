 // ННГУ, ИИТММ, Курс "Алгоритмы и структуры данных"
//
// tcomputecentermodel.cpp
//
// Copyright (c) Пинежанин Е.С.

#include "tpolynom.h"
#include <sstream>
#include <cctype>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <regex>


TPolynom operator*(double coeff, const TPolynom& polynom) {
    return polynom * coeff;   
}

bool TPolynom::operator==(const TPolynom& polynom) const {
    if (data.size() != polynom.data.size()) return false;

    auto it1 = data.begin();
    auto it2 = polynom.data.begin();

    while (it1 != data.end()) {
        if (it1->powers != it2->powers) return false;
        if (std::abs(it1->coeff - it2->coeff) > 1e-9) return false;
        ++it1; ++it2;
    }
    return true;
}

bool TPolynom::operator!=(const TPolynom& polynom) const {
    return !(*this == polynom);
}


// Вспомогательные функции для работы со степенями
namespace {
    /**
     * Степени переменных кодируются в одно число по формуле:
     * powers = a*100 + b*10 + c
     * где a - степень x, b - степень y, c - степень z
     */
    const int X_MULTIPLIER = 100;
    const int Y_MULTIPLIER = 10;
    const int Z_MULTIPLIER = 1;

    //Создание свернутого индекса из степеней
    int MakePowers(int x_deg, int y_deg, int z_deg) {
        return x_deg * X_MULTIPLIER + y_deg * Y_MULTIPLIER + z_deg;
    }
    //Извлечение степеней из свернутого индекса
    void ExtractPowers(int powers, int& x_deg, int& y_deg, int& z_deg) {
        x_deg = powers / X_MULTIPLIER;
        y_deg = (powers % X_MULTIPLIER) / Y_MULTIPLIER;
        z_deg = powers % Y_MULTIPLIER;
    }

    //Проверка, является ли символ допустимой переменной
    bool IsVariable(char ch) {
        return ch == 'x' || ch == 'y' || ch == 'z';
    }


    //Парсинг числа из строки  <значение, новая_позиция>

    std::pair<double, size_t> ParseNumber(const std::string& str, size_t pos) {
        double result = 0;
        bool negative = false;
        bool has_decimal = false;
        double decimal_place = 0.1;

        if (pos < str.length() && str[pos] == '-') {
            negative = true;
            pos++;
        }
        else if (pos < str.length() && str[pos] == '+') {
            pos++;
        }

        while (pos < str.length() && (std::isdigit(str[pos]) || str[pos] == '.')) {
            if (str[pos] == '.') {
                has_decimal = true;
                pos++;
                continue;
            }

            if (!has_decimal) {
                result = result * 10 + (str[pos] - '0');
            }
            else {
                result = result + (str[pos] - '0') * decimal_place;
                decimal_place /= 10;
            }
            pos++;
        }

        // Если после знака нет числа, значит коэффициент = 1
        if (result == 0 && pos < str.length() && IsVariable(str[pos])) {
            result = 1;
        }

        return { negative ? -result : result, pos };
    }

    //Парсинг степени переменной pair<степень, новая_позиция>
    std::pair<int, size_t> ParseDegree(const std::string& str, size_t pos) {
        int degree = 1; // По умолчанию степень 1

        if (pos < str.length() && str[pos] == '^') {
            pos++;
            degree = 0;
            while (pos < str.length() && std::isdigit(str[pos])) {
                degree = degree * 10 + (str[pos] - '0');
                pos++;
            }
        }

        return { degree, pos };
    }

    //Сравнение мономов для сортировки (по убыванию суммарной степени)
    bool CompareTerms(const TTerm& a, const TTerm& b) {
        int a_x, a_y, a_z, b_x, b_y, b_z;
        ExtractPowers(a.powers, a_x, a_y, a_z);
        ExtractPowers(b.powers, b_x, b_y, b_z);

        // Сортируем по убыванию суммарной степени
        int a_sum = a_x + a_y + a_z;
        int b_sum = b_x + b_y + b_z;

        if (a_sum != b_sum) {
            return a_sum > b_sum;
        }

        // При равной сумме - лексикографически по x, y, z
        if (a_x != b_x) return a_x > b_x;
        if (a_y != b_y) return a_y > b_y;
        return a_z > b_z;
    }

    //Приведение подобных членов и сортировка
    void SimplifyPolynom(std::list<TTerm>& data) {
        if (data.empty()) return;

        // Сортируем для группировки
        data.sort(CompareTerms);

        // Объединяем подобные члены
        auto it = data.begin();
        while (it != data.end()) {
            auto next = std::next(it);
            if (next != data.end() && it->powers == next->powers) {
                it->coeff += next->coeff;
                data.erase(next);
                // Проверяем, не обнулился ли коэффициент
                if (std::abs(it->coeff) < 1e-10) {
                    it = data.erase(it);
                }
            }
            else {
                ++it;
            }
        }
    }
}

 // Конструкторы
 

TPolynom::TPolynom() {
    // data инициализируется автоматически
}
 
void TPolynom::ParsePolynom(const std::string& expr) {
    size_t pos = 0;
    bool first_term = true;

    while (pos < expr.length()) {
        double coeff = 1.0;
        int x_deg = 0, y_deg = 0, z_deg = 0;

        // Обработка знака
        if (!first_term) {
            if (expr[pos] == '+') {
                coeff = 1.0;
                pos++;
            }
            else if (expr[pos] == '-') {
                coeff = -1.0;
                pos++;
            }
        }
        else {
            if (expr[pos] == '-') {
                coeff = -1.0;
                pos++;
            }
            first_term = false;
        }

        // Парсим коэффициент
        if (pos < expr.length() && (std::isdigit(expr[pos]) || expr[pos] == '.')) {
            double num;
            std::tie(num, pos) = ParseNumber(expr, pos);
            coeff *= num;
        }

        // Парсим переменные
        while (pos < expr.length() && IsVariable(expr[pos])) {
            char var = expr[pos];
            pos++;

            int degree;
            std::tie(degree, pos) = ParseDegree(expr, pos);

            switch (var) {
            case 'x': x_deg = degree; break;
            case 'y': y_deg = degree; break;
            case 'z': z_deg = degree; break;
            }
        }

        // Добавляем моном
        if (std::abs(coeff) > 1e-10) {
            TTerm term;
            term.coeff = coeff;
            term.powers = MakePowers(x_deg, y_deg, z_deg);
            data.push_back(term);
        }
    }
}


 
bool TPolynom::IsCorrect() const {
    for (const auto& term : data) {
        int x_deg, y_deg, z_deg;
        ExtractPowers(term.powers, x_deg, y_deg, z_deg);

        // Проверка степеней (0-9 по условию лабораторной работы)
        if (x_deg < 0 || x_deg > 9 ||
            y_deg < 0 || y_deg > 9 ||
            z_deg < 0 || z_deg > 9) {
            return false;
        }

        // Проверка коэффициента (не бесконечность и не NaN)
        if (std::isinf(term.coeff) || std::isnan(term.coeff)) {
            return false;
        }
    }
    return true;
}
 // Установка полинома из строки
void TPolynom::SetPolynom(const std::string& polynom) {
    data.clear();

    if (polynom.empty()) {
        return;
    }

    std::string expr = polynom;
    expr.erase(std::remove_if(expr.begin(), expr.end(), ::isspace), expr.end());

    if (expr.empty()) {
        return;
    }

    ParsePolynom(expr);

    SimplifyPolynom(data);
}
 
// Арифметические операции
TPolynom TPolynom::operator+(const TPolynom& polynom) const {
    TPolynom result;
    result.data = data;

    for (const auto& term : polynom.data) {
        result.data.push_back(term);
    }

    SimplifyPolynom(result.data);
    return result;
}

 TPolynom TPolynom::operator-(const TPolynom& polynom) const {
     TPolynom result;
     result.data = data;

     for (const auto& term : polynom.data) {
         TTerm neg_term;
         neg_term.coeff = -term.coeff;
         neg_term.powers = term.powers;
         result.data.push_back(neg_term);
     }

     SimplifyPolynom(result.data);
     return result;
 }

 TPolynom TPolynom::operator*(const TPolynom& polynom) const {
     TPolynom result;

     for (const auto& term1 : data) {
         for (const auto& term2 : polynom.data) {
             TTerm product;
             product.coeff = term1.coeff * term2.coeff;

             int x1, y1, z1, x2, y2, z2;
             ExtractPowers(term1.powers, x1, y1, z1);
             ExtractPowers(term2.powers, x2, y2, z2);

             int new_x = x1 + x2;
             int new_y = y1 + y2;
             int new_z = z1 + z2;

             if (new_x <= 9 && new_y <= 9 && new_z <= 9) {
                 product.powers = MakePowers(new_x, new_y, new_z);
                 if (std::abs(product.coeff) > 1e-10) {
                     result.data.push_back(product);
                 }
             }
         }
     }

     SimplifyPolynom(result.data);
     return result;
 }

 TPolynom TPolynom::operator*(double coeff) const {
     TPolynom result;

     if (std::abs(coeff) < 1e-10) {
         return result;
     }

     for (const auto& term : data) {
         TTerm new_term;
         new_term.coeff = term.coeff * coeff;
         new_term.powers = term.powers;

         if (std::abs(new_term.coeff) > 1e-10) {
             result.data.push_back(new_term);
         }
     }

     SimplifyPolynom(result.data);
     return result;
 }
void TPolynom::Add(const std::string& monom)   {
    TPolynom temp;
    temp.SetPolynom(monom);

    for (const auto& term : temp.data) {
        data.push_back(term);
    }

    SimplifyPolynom(data);
}

 
void TPolynom::Delete(size_t pos) {
    if (pos < data.size()) {
        auto it = data.begin();
        std::advance(it, pos);
        data.erase(it);
    }
}
 
double TPolynom::Calculate(double x, double y, double z) const {
    double result = 0.0;

    for (const auto& term : data) {
        int x_deg, y_deg, z_deg;
        ExtractPowers(term.powers, x_deg, y_deg, z_deg);

        double term_value = term.coeff *
            std::pow(x, x_deg) *
            std::pow(y, y_deg) *
            std::pow(z, z_deg);
        result += term_value;
    }

    return result;
}
 
std::istream& operator>>(std::istream& istr, TPolynom& polynom) {
    std::string input;
    std::getline(istr, input);
    polynom.SetPolynom(input);
    return istr;
}

 
std::ostream& operator<<(std::ostream& ostr, const TPolynom& polynom) {
    if (polynom.data.empty()) {
        ostr << "0";
        return ostr;
    }

    bool first = true;

    for (const auto& term : polynom.data) {
        int x_deg, y_deg, z_deg;
        ExtractPowers(term.powers, x_deg, y_deg, z_deg);

        double coeff = term.coeff;

        // Печать знака
        if (!first) {
            if (coeff > 0) {
                ostr << " + ";
            }
            else {
                ostr << " - ";
            }
        }
        else {
            if (coeff < 0) {
                ostr << "-";
            }
            first = false;
        }

        double abs_coeff = std::abs(coeff);

        // Печать коэффициента (опускаем 1, если есть переменные)
        bool has_variables = (x_deg > 0 || y_deg > 0 || z_deg > 0);
        if (abs_coeff != 1.0 || !has_variables) {
            // Если коэффициент целый, печатаем без десятичной точки
            if (abs_coeff == static_cast<int>(abs_coeff)) {
                ostr << static_cast<int>(abs_coeff);
            }
            else {
                ostr << abs_coeff;
            }
        }

        // Печать переменных и степеней
        if (x_deg > 0) {
            ostr << "x";
            if (x_deg > 1) ostr << "^" << x_deg;
        }
        if (y_deg > 0) {
            ostr << "y";
            if (y_deg > 1) ostr << "^" << y_deg;
        }
        if (z_deg > 0) {
            ostr << "z";
            if (z_deg > 1) ostr << "^" << z_deg;
        }
    }

    return ostr;
}