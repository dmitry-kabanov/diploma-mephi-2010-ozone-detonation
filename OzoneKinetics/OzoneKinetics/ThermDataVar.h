#ifndef THERM_DATA_VAR_H
#define THERM_DATA_VAR_H

class ThermDataVar  
{
public:
	/**
     * Выделяет память в куче под массивы термоданных.
     *
     * @return Всегда возвращает ноль.
     */
    int AllocateMemoryForTemperatureRange();
    /**
     * Конструктор класса.
     */
	ThermDataVar();
    /**
     * Деструктор класса.
     */
	virtual ~ThermDataVar();
	/**
     * Количество температурных интервалов.
     */
    int n;
    /**
     * Нижний предел температуры.
     */
	double *Tlow;
    /**
     * Верхний предел температуры.
     */
	double *Tup;
    /**
     * Набор коэффициентов для каждого вещества.
     */
	double **a;
};

#endif
