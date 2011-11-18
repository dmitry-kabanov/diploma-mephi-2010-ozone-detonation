#include <iostream>
#include "AirlineTicketLib/AirlineTicket.h"

using namespace std;
using namespace AirlineTicketLib;

int main()
{
	AirlineTicket myTicket;
	myTicket.setPassengerName("Sherman T. Socketwrench");
	myTicket.setNumberOfMiles(700);
	int cost = myTicket.calculatePriceInDollars();
	cout << "This ticket costs $" << cost << endl;

	AirlineTicket *myTicket2 = new AirlineTicket();
	myTicket2->setPassengerName("Laudimore M. Hallidue");
	myTicket2->setNumberOfMiles(2000);
	myTicket2->setHasEliteSuperRewardsStatus(true);
	int cost2 = myTicket2->calculatePriceInDollars();
	cout << "And this ticket costs $" << cost2 << endl;
	delete myTicket2;

	return EXIT_SUCCESS;
}