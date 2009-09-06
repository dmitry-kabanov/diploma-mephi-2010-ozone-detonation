#include "Elements.h"

Elements::Elements()
{
	ElementsList[0] = 'H';
	ElementsList[1] = 'O';
	ElementsList[2] = 'C';
	ElementsList[3] = 'N';

	N = 4;
	int i;
	for (i = 0; i < N; i++)
		nBeta[i] = 0;
}

Elements::~Elements()
{

}
