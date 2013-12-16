#include "Breaker.h"

void Breaker::switchIt()
{
	if (state == 1)
	{
		state = 0;
	} 
	else
	{
		state = 1;
	}
}