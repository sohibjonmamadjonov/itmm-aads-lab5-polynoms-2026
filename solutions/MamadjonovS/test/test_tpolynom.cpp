
#include "tpolynom.h"

#include <gtest.h>
 
TEST(TPolynom, SimpleAddition)
{
    TPolynom P;
    P.SetPolynom("3x^5y^2z^5 - 5x^4y^3z^3 + 7x^3y^5z");

    TPolynom Q;
    Q.SetPolynom("4x^3y^2z^6 - 6x^2yz^8");

    TPolynom expected;
    expected.SetPolynom("3x^5y^2z^5 - 5x^4y^3z^3 + 7x^3y^5z + 4x^3y^2z^6 - 6x^2yz^8");

    TPolynom result = P + Q;

    EXPECT_TRUE(result == expected)
        << "Ошибка в тесте 1: Неправильное сложение полиномов\n"
        << "P = 3x^5y^2z^5 - 5x^4y^3z^3 + 7x^3y^5z\n"
        << "Q = 4x^3y^2z^6 - 6x^2yz^8";
     
}
 
TEST(TPolynom, AdditionWithDifferentOrder)
{
    TPolynom P;
    P.SetPolynom("3x^5y^2z^5 - 5x^4y^3z^3 + 7x^3y^5z");

    TPolynom Q;
    Q.SetPolynom("4x^7y^2z^6 - 6x^6yz^8");

    TPolynom expected;
    expected.SetPolynom("4x^7y^2z^6 - 6x^6yz^8 + 3x^5y^2z^5 - 5x^4y^3z^3 + 7x^3y^5z");

    TPolynom result = P + Q;


    EXPECT_TRUE(result == expected)
        << "Ошибка в тесте 2: Неправильная сортировка при сложении\n"
        << "Члены с большими степенями должны быть в начале";
     
}
 
TEST(TPolynom, LikeTermsCombination)
{
    TPolynom P;
    P.SetPolynom("3x^5y^2z^5 - 5x^4y^3z^3 + 7x^3y^5z");

    TPolynom Q;
    Q.SetPolynom("4x^5y^2z^5 + 5x^4y^3z^3");

    TPolynom expected;
    expected.SetPolynom("7x^5y^2z^5 + 7x^3y^5z");

    TPolynom result = P + Q;

    EXPECT_TRUE(result == expected)
        << "Ошибка в тесте 3: Неправильное приведение подобных членов\n"
        << "Коэффициенты должны складываться: 3+4=7, -5+5=0 (член удаляется)";
 
}

 
TEST(TPolynom, ComplexOrdering)
{
    TPolynom P;
    P.SetPolynom("3x^5y^2z^5 - 5x^4y^3z^3 + 7x^7y^5z");

    TPolynom Q;
    Q.SetPolynom("4x^6y^2z^6 - 6x^2yz^8");

    TPolynom expected;
    expected.SetPolynom("7x^7y^5z + 4x^6y^2z^6 + 3x^5y^2z^5 - 5x^4y^3z^3 - 6x^2yz^8");

    TPolynom result = P + Q;

    EXPECT_TRUE(result == expected)
        << "Ошибка в тесте 4: Неправильная сортировка\n"
        << "Члены должны быть отсортированы по убыванию суммарной степени";
 
}

 
TEST(TPolynom, CompleteCancellation)
{
    TPolynom P;
    P.SetPolynom("3x^5y^2z^5 - 5x^4y^3z^3 + 7x^7y^5z");

    TPolynom Q;
    Q.SetPolynom("-3x^5y^2z^5 + 5x^4y^3z^3 - 7x^7y^5z");

    TPolynom expected;
    expected.SetPolynom("0");

    TPolynom result = P + Q;
 
    EXPECT_TRUE(result == expected)
        << "Ошибка в тесте 5: При сложении противоположных полиномов\n"
        << "Должен получаться нулевой полином '0'";
}

 
TEST(TPolynom, ConstructorAndIO)
{
    TPolynom p1;
    std::stringstream ss1;
    ss1 << p1;
    EXPECT_EQ(ss1.str(), "0") << "Пустой полином должен выводиться как '0'";

    TPolynom p2;
    p2.SetPolynom("2x^2 + 3y^2 - 4z^2");

    EXPECT_TRUE(p2.IsCorrect()) << "Полином 2x^2 + 3y^2 - 4z^2 должен быть корректным";
}
 
TEST(TPolynom, Multiplication)
{
    TPolynom P;
    P.SetPolynom("2x + 3y");

    TPolynom Q;
    Q.SetPolynom("4x + 5y");

    TPolynom expected;   
    expected.SetPolynom("8x^2 + 22xy + 15y^2");

    TPolynom result = P * Q;

    EXPECT_TRUE(result == expected) 
        << "Ошибка в тесте умножения: (2x+3y)*(4x+5y) должно быть 8x^2+22xy+15y^2";
}
 
TEST(TPolynom, Calculation)
{
    TPolynom P;
    P.SetPolynom("2x + 3y + 4z");

    double result = P.Calculate(1.0, 2.0, 3.0);

    // 2*1 + 3*2 + 4*3 = 2 + 6 + 12 = 20
    EXPECT_DOUBLE_EQ(result, 20.0)
        << "Вычисление значения полинома 2x+3y+4z в точке (1,2,3) должно дать 20";
}

 
TEST(TPolynom, AddMonom)
{
    TPolynom P;
    P.SetPolynom("2x^2 + 3y^2");

    P.Add("4z^2");

    TPolynom expected;
    expected.SetPolynom("2x^2 + 3y^2 + 4z^2");

    EXPECT_TRUE(P == expected)  
        << "После добавления 4z^2 должен получиться полином 2x^2 + 3y^2 + 4z^2";
}
 
TEST(TPolynom, DeleteMonom)
{
    TPolynom P;
    P.SetPolynom("2x^2 + 3y^2 + 4z^2");

    P.Delete(1);
 
    TPolynom expected;
    expected.SetPolynom("2x^2 + 4z^2");

    

    EXPECT_TRUE(P == expected)  
        << "Должен получиться полином '2x^2 + 4z^2'";
}

TEST(TPolynom, InequalityOperator) {
    TPolynom p1, p2;
    p1.SetPolynom("6x^4 + 2y^4 + 5z^3");
    p2.SetPolynom("6x^4 + 2y^4 + 2z^3");

    EXPECT_TRUE(p1 != p2) << "Разные полиномы должны быть неравны";
}

TEST(TPolynom, EqualityDifferentOrder) {
    TPolynom p1, p2;
    p1.SetPolynom("6x^4 + 2y^4 + 5z^3");
    p2.SetPolynom("2y^4 + 5z^3 + 6x^4");   

    EXPECT_TRUE(p1 == p2) << "Полиномы должны быть равны независимо от порядка членов";
}


TEST(TPolynom, DoubleCoefficientsAddition)
{
    TPolynom P;
    P.SetPolynom("3.5x^2 + 2.7y^2");

    TPolynom Q;
    Q.SetPolynom("1.3x^2 + 4.2y^2");

    TPolynom expected;
    expected.SetPolynom("4.8x^2 + 6.9y^2");

    TPolynom result = P + Q;

    EXPECT_TRUE(result == expected)
        << "3.5x^2 + 2.7y^2 + 1.3x^2 + 4.2y^2 = 4.8x^2 + 6.9y^2";
}

TEST(TPolynom, DoubleCoefficientsMultiplication)
{
    TPolynom P;
    P.SetPolynom("2.5x + 1.5y");

    TPolynom Q;
    Q.SetPolynom("3.0x + 2.0y");

    TPolynom expected;
    expected.SetPolynom("7.5x^2 + 9.5xy + 3y^2");   

    TPolynom result = P * Q;

    EXPECT_TRUE(result == expected)
        << "(2.5x+1.5y)*(3.0x+2.0y) = 7.5x^2 + 9.5xy + 3y^2";
}

TEST(TPolynom, DoubleNumberMultiplication)
{
    TPolynom P;
    P.SetPolynom("2x^2 + 3y^2");

    TPolynom result = 5.34 * P;

    TPolynom expected;
    expected.SetPolynom("10.68x^2 + 16.02y^2");

    EXPECT_TRUE(result == expected)
        << "5.34 * (2x^2 + 3y^2) = 10.68x^2 + 16.02y^2";
}

TEST(TPolynom, DoubleCoefficientsCalculation)
{
    TPolynom P;
    P.SetPolynom("2.5x + 1.5y + 0.5z");

    double result = P.Calculate(2.0, 3.0, 4.0);

    EXPECT_DOUBLE_EQ(result, 11.5)
        << "2.5x+1.5y+0.5z в точке (2,3,4) = 2.5*2 + 1.5*3 + 0.5*4 = 5 + 4.5 + 2 = 11.5";
}