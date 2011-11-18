#ifndef AIRLINETICKET_AIRLINETICKETLIB_AIRLINETICKET_H
#define AIRLINETICKET_AIRLINETICKETLIB_AIRLINETICKET_H

#include <string>

namespace AirlineTicketLib
{

	class AirlineTicket
	{
	public:
		AirlineTicket();
		~AirlineTicket();

		int calculatePriceInDollars();

		std::string getPassengerName();
		void setPassengerName(std::string name);
		int getNumberOfMiles();
		void setNumberOfMiles(int miles);
		bool getHasEliteSuperRewardsStatus();
		void setHasEliteSuperRewardsStatus(bool status);

	private:
		std::string passengerName;
		int numberOfMiles;
		bool hasEliteSuperRewardsStatus;
	};

}

#endif