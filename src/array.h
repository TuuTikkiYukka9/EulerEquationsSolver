#pragma once
#include "Solver.h"

template <typename T>
class Array
{
private:
	int _size;
	T *arr;
public:
	Array(int size) {
		_size = size;
		arr = ((_size > 0) ? new T[_size] : nullptr);
	};
	~Array() {
		delete[] arr;
	};

	T& operator[](int index) {
		if (index >= 0 && index < _size) {
			return arr[index];
		}
		else {
			throw "Index outside the bounds of the array!";
		}
	}

	int length() {
		return _size;
	};


	void setFirst(T value) {
		arr[0] = value;
	};

	void setLast(T value) {
		arr[_size - 1] = value;
	}


	T max() {
		T maxValue = arr[0];
		for (int i = 1; i < _size; i++) {
			if (arr[i] > maxValue) maxValue = arr[i];
		}
		return maxValue;
	}
};

