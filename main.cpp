#include "mainwindow.h"
#include "wrapper.h"
#include <QApplication>
#include <iostream>
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    std::vector<double> x = wrapper(1);
    for(auto i: x)
        std::cout << i << std::endl;
    return a.exec();
}
