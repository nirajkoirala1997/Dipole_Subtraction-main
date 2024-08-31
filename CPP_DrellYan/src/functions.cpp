// rectangle.cpp
#include "rectangle.h"
#include <iostream>

Rectangle::Rectangle(int w, int h) : width(w), height(h) {}

int Rectangle::area() {
    return width * height;
}

void Rectangle::display() {
    std::cout << "Width: " << width << ", Height: " << height << std::endl;
}

void Rectangle::setWidth(int w) {
    width = w;
}

void Rectangle::setHeight(int h) {
    height = h;
}

int Rectangle::getWidth() {
    return width;
}

int Rectangle::getHeight() {
    return height;
}



//
//
//
//
//#include <iostream>
//#include "rectangle.h"
//using namespace std;
//
//class Rectangle {
//private:
//    // Private data members (only accessible within the class)
//    int width;
//    int height;
//
//public:
//    // Constructor
//    Rectangle(int w, int h) : width(w), height(h) {}
//
//    // Member function to calculate area
//    int area() {
//        return width * height;
//    }
//
//    // Member function to display the dimensions
//    void display() {
//        cout << "Width: " << width << ", Height: " << height << endl;
//    }
//
//    // Setter for width
//    void setWidth(int w) {
//        width = w;
//    }
//
//    // Setter for height
//    void setHeight(int h) {
//        height = h;
//    }
//
//    // Getter for width
//    int getWidth() {
//        return width;
//    }
//
//    // Getter for height
//    int getHeight() {
//        return height;
//    }
//};
