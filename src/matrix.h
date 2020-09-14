#pragma once
#include "array.h"

template <typename T>
class Matrix
{
private:
	int _height;
	int _width;
	Array<Array<T>*>* mat;
public:
	Matrix(int height, int width) {
		_height = height; 
		_width = width;
		mat = new Array<Array<T>*>(_height);
		for (int i = 0; i < mat->length(); i++)
		{
			(*mat)[i] = new Array<T>(_width);
		}
	};

	Array<T>& operator[](int index) {
		if (index >= 0 && index < _height) {
			return (*(*mat)[index]);
		}
		else {
			throw "Index outside the bounds of the array!";
		}
	};

	void SetAllElements(T value) {
		for (int i = 0; i < mat->length(); i++)
		{
			for (int j = 0; j < ((*mat)[i])->length(); j++)
			{
				(*(*mat)[i])[j] = value;
			}
		}
	}

	~Matrix() {
		for (int i = 0; i < mat->length(); i++)
		{
			delete (*mat)[i];
		}
		delete mat;
	};

	int height() {
		return _height;
	}

	int width() {
		return _width;
	}
	
};

