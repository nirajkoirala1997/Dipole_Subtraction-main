#include <iostream>
#include "rectangle.h"  // Include the Rectangle class declaration

using namespace std;

int main() {
    Rectangle rect(5, 10);
    rect.display();
    cout << "Area: " << rect.area() << endl;

    rect.setWidth(7);
    rect.setHeight(12);

    rect.display();
    cout << "New Area: " << rect.area() << endl;

    return 0;
}

