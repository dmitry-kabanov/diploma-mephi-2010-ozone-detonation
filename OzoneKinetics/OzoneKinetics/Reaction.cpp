/**
 * @file
 *
 * @author  Кабанов Дмитрий <kabanovdmitry@gmail.com
 * @version %I%
 *
 * @section DESCRIPTION
 *
 * Реализация класса Reaction.
 */
#include "Reaction.h"

Reaction::Reaction()
{
}

Reaction::~Reaction()
{
    delete [] reagents;
    delete [] products;
}