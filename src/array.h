#pragma once
#include "Solver.h"

template <typename T>
class Array
{
private:
	int _size;
	T *arr;
public:
	Array() {
		_size = 0;
		arr = ((_size > 0) ? new T[_size] : nullptr);
	};

	Array(int size) {
		_size = size;
		arr = ((_size > 0) ? new T[_size] : nullptr);
	};
	~Array() {
		if (_size > 0) {
			delete[] arr;
		}
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

	//void resize(int size) {
	//	T* newArr = new T[size];
	//	memcpy(newArr, arr, (size > _size ? _size : size) * sizeof(T));
	//	_size = size;
	//	delete[] arr;
	//	arr = newArr;
	//}

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

