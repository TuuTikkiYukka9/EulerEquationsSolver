#pragma once
#include <cmath>
#include <string>
#include <iostream>

double to_double(std::string str) {
	int znak = 1;
	int razrad = 10;

	if (str.size() == 0) {
		return 0.0;
	}
	else if (str[0] == '-') {
		znak = -1;
		std::string buf = "";
		for (int i = 1;i < str.size();i++) {
			buf += str[i];
		}
		str = buf;
		//buf = "";
		//std::cout << "minus : " << str << "\n";
	}

	std::string buf = "";
		for (int i = 0;i < str.size();i++) {
			if (str[i] == '.') razrad = pow(10, (str.size() - i));
			else buf += str[i];
		}
		str = buf;
		//std::cout << "del '.' : " << str << "\n";


	double res = 0;
	int len = str.size();
	for (int i = 0;i < len;i++) {
		if (str[i] == '1') res += 1 * (pow(10, (len-i)));
		else if (str[i] == '2') res += 2 * (pow(10, (len-i)));
		else if (str[i] == '3') res += 3 * (pow(10, (len-i)));
		else if (str[i] == '4') res += 4 * (pow(10, (len-i)));
		else if (str[i] == '5') res += 5 * (pow(10, (len-i)));
		else if (str[i] == '6') res += 6 * (pow(10, (len-i)));
		else if (str[i] == '7') res += 7 * (pow(10, (len-i)));
		else if (str[i] == '8') res += 8 * (pow(10, (len-i)));
		else if (str[i] == '9') res += 9 * (pow(10, (len-i)));
	}
	return (res / razrad)*znak;
}
