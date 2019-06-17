#include <iostream>
#include <string.h>

#include "gps-sdr-sim/gpssim.h"
#include "zmq.hpp"
#define MAX_LENGTH 10


int main() {
	struct structB{
	    int a;
	    char c[MAX_LENGTH];
	};

	//struct structB messageB;
	//messageB.a=0;
	//strcpy(messageB.c,"aa"); // #include<cstring> or use std::copy from <algorithm>
    gpstime_t g1;
    g1.week = 1989;//weeks from 1980
    g1.sec = 172800;//2-12-2018


	const int length = sizeof(int) + sizeof(double) + 1;
	zmq::message_t msg (length);
	memcpy (msg.data(), &g1, length);
    std::cout << "Message sent" << std::endl;
	//socket->send(msg);

	return 0;
}
