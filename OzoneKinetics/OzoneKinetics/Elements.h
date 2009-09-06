#ifndef ELEMENTS_H
#define ELEMENTS_H

class Elements  
{
public:
    /**
     * Конструктор класса.
     */
	Elements();
    /**
     * Деструктор класса.
     */
	virtual ~Elements();
    /**
     * Массив химических элементов (C, O, H, N).
     */
	char ElementsList[100];
    /**
     * Что-то непонятное.
     */
	int nBeta[100];
    /**
     * Количество химических элементов.
     */
	int N;
};

#endif
