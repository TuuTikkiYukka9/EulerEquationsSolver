#pragma once
#include "solver.h"

template <typename T>
class Array
{
private:
	int _size;
	T *arr;
public:
	Array() : Array(0) {};

	Array(int size) {
		_size = size;
		arr = ((_size > 0) ? new T[_size] : nullptr);
	};
	Array(int size, T defualtValue) : Array(size) {
		for (int i = 0; i < _size; i++) {
			arr[i] = defualtValue;
		}
	};
	~Array() {
		if (_size > 0 && arr != nullptr) {
			delete[] arr;
		}
		arr = nullptr;
		_size = 0;
	};

	T& operator[](int index) {
		if (index >= 0 && index < _size) {
			return arr[index];
		}
		else {
			throw "Index outside the bounds of the array!";
		}
	}

	T& item(int index) {
		if (index >= 0 && index < _size) {
			return arr[index];
		}
		else {
			throw "Index outside the bounds of the array!";
		}
	}

	Array& operator=(const Array& right) {
		//проверка на самоприсваивание
		if (this == &right) {
			return *this;
		}
		delete[] arr;
		_size = right._size;
		arr = ((_size > 0) ? new T[_size] : nullptr);
		for (int i = 0; i < _size; i++){
			arr[i] = right.arr[i];
		}
		return *this;
	}

	void copyValues(const Array& right) {
		if (_size != right._size) {
			throw "Arrays must have the same length!";
		}
		for (int i = 0; i < _size; i++) {
			arr[i] = right.arr[i];
		}
	}

	void print() {
		for (int i = 0; i < _size; i++) {
			std::cout << arr[i] << ", ";
		}
		std::cout << "\n";
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

