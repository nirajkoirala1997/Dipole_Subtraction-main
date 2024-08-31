// rectangle.h
#ifndef RECTANGLE_H
#define RECTANGLE_H

class Rectangle {
private:
    int width;
    int height;

public:
    Rectangle(int w, int h);
    int area();
    void display();
    void setWidth(int w);
    void setHeight(int h);
    int getWidth();
    int getHeight();
};

#endif // RECTANGLE_H

