#ifndef MIXTURE_H
#define MIXTURE_H

#include <string>
#include "ThermDataVar.h"
#include "ReactionDataVar.h"
#include "Elements.h"

#define Rgas 8.31434
#define Rgas_ext 8314.34
#define Na 6.022045E23

class Mixture
{
public:
    char ErrorBuffer[250];
    /**
     * Читает файл с веществами.
     * 
     * @param fileName имя файла с веществами
     */
	int readFileOfSubstances(const char *fileName);
    /**
     * Читает файл с реакциями.
     *
     * @param fileName имя файла с реакциями
     */
	int readFileOfReactions(const char *fileName);
    /**
     * Конструктор класса.
     */
	Mixture();
    /**
     * Деструктор класса.
     */
	virtual ~Mixture();
    /**
     * Количество веществ в смеси.
     */
	int nSpecies;
    /**
     * Количество реакций.
     */
	int nReactions;
    /**
     * Термодинамические данные для каждого из веществ.
     */
	ThermDataVar *thermDatVar;
	double *Gibbs_Ef;
	double *MU;
	double *H0;
    /**
     * Названия веществ.
     */
    std::string *namesOfSubstances;
	Elements *Atom;
	ReactionDataVar *reactionsDatVar;
    /**
     * Имя файла ошибок.
     */
	char *nameOfFileOfErrors;
};

#endif