#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <complex>
#include <Eigen/Dense>
QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void makePlot();
    void on_Simulate_clicked();

private:
    Ui::MainWindow *ui;
    std::vector<double> simulate();
    void run();
};
#endif // MAINWINDOW_H
