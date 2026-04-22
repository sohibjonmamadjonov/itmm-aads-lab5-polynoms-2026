// ННГУ, ИИТММ, Курс "Алгоритмы и структуры данных"
//
// polynom.h
//
// Copyright (c) Пинежанин Е.С.

#ifndef __TPolynom_H__
#define __TPolynom_H__

#include <list>
#include <string>

struct TTerm
{
  double coeff;
  int powers;
};

class TPolynom {
private:
	std::list<TTerm> data;

	void ParsePolynom(const std::string& expr);

public:
	TPolynom();

	bool operator==(const TPolynom& other) const;
	bool operator!=(const TPolynom& other) const;

	bool IsCorrect() const; // возвращает коррекность полинома (false - если на вход подали некорректное выражение)

	void SetPolynom(const std::string& polynom);

	TPolynom operator+(const TPolynom& polynom) const;
	TPolynom operator-(const TPolynom& polynom) const;
	TPolynom operator*(const TPolynom& polynom) const;
	TPolynom operator*(double coeff) const;

	void Add(const std::string& monom);  
	void Delete(size_t pos);  

	double Calculate(double x, double y, double z) const;

	friend std::istream& operator>>(std::istream& istr, TPolynom& polynom);
	friend std::ostream& operator<<(std::ostream& ostr, const TPolynom& polynom);
};

TPolynom operator*(double coeff, const TPolynom& polynom);

#endif
