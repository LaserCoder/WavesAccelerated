#include "mainwindow.h"
#include "qcustomplot.h"
#include "./ui_mainwindow.h"
#include <cmath>
#include <iostream>
#include <QTimer>
#include <Eigen/Dense>
#include <complex>
#include <chrono>
#include <thread>
#include <fftw3.h>

enum lineshape {Lorentzian, Parabolic};
enum phi {Fundamental, HarmonicN, CosineN, Rand};
struct Params{
    double J;
    double kpp;
    double gammaK;
    double R1;
    double R2;
    double lm0;
    double Lc;
    double df;

    double mu;
    double T1;
    double T2;
    double n;
    double Lmod;

    double aw;
    double Gamma;

    double h;
    int dtperTr;
    int numTr;
    double Crnt;
    bool gc;
    double gcfac;
    lineshape ls;
    bool useN2N3;
    phi initphi;

    bool plotprogress;
    bool plottheory;
    bool shifttomin;
};

int t = 1;
Params sim_params;
auto out = Eigen::Array<std::complex<double>, 1, 5001>(0);

void refreshGraph(Ui::MainWindow* ui)
{
    Eigen::ArrayXf x = Eigen::ArrayXf::LinSpaced(5001, 0, 1);
    Eigen::Array<double, 5001, 1> y = out.real();
    auto a = QVector<double>(x.data(),x.data() + x.size());
    auto b = QVector<double>(y.data(),y.data() + y.size());
    ui->rightPlot->graph(0)->setData(a,b);
    ui->rightPlot->replot();
    ui->rightPlot->update();
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    MainWindow::makePlot();
    
    
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::makePlot()
{
    ui->rightPlot->addGraph();
    ui->rightPlot->xAxis->setLabel("x");
    ui->rightPlot->yAxis->setLabel("y");
    ui->rightPlot->xAxis->setRange(0, 1);
    ui->rightPlot->yAxis->setRange(-6e6,6e6 );

    ui->leftPlot->addGraph();
    ui->leftPlot->xAxis->setLabel("x");
    ui->leftPlot->yAxis->setLabel("y");
    ui->leftPlot->xAxis->setRange(0, 1);
    ui->leftPlot->yAxis->setRange(-6e6,6e6 );
    std::cout << "Hello World!" << std::endl;
}

std::complex<double> ii (0, 1.0);
using Eigen::Array;
void MainWindow::simulate()
{
      double J = sim_params.J;
//    double R1 = sim_params.R1;
//    double R2 = sim_params.R2;
//    double c = 299792458;
//    double n = sim_params.n;
//    double kpp = sim_params.kpp;
//    double aw = sim_params.aw;
//    double gammaK = sim_params.gammaK;
//    double g0 = 819.429;
//    double gc = sim_params.gc;
//    double T1 = sim_params.T1;
//    double T2 = sim_params.T2;
//    double dz =  8.0000e-07;
//    double dt = 8.8061e-15;
//    double Psat = 8.189848267404146e+12;
    auto start = std::chrono::high_resolution_clock::now();
    double R1 = 0.09;
    double R2 = 1;
    double c = 299792458;
    double n = 3.3;
    double kpp = -2.000000000000000e-24;
    double aw = 400;
    double gammaK = 0;
    double g0 = 819.4299;
    double gc = 1;
    double T1 = 4.0000e-13;
    double T2 = 5.0000e-14;
    double dz =  8.0000e-07;
    double dt = 8.8061e-15;
    double Psat = 8.189848267404146e+12;
    Array<std::complex<double>, 1 , 5001> Ep;
    Array<std::complex<double>, 1 , 5001> Em;
    Array<std::complex<double>, 1 , 5001> Pp;
    Array<std::complex<double>, 1 , 5001> Pm;
    Array<std::complex<double>, 1 , 5003> tmpdz2;
    Array<std::complex<double>, 1 , 5002> tmpdz;
    Array<std::complex<double>, 1 , 5001> dEpdz;
    Array<std::complex<double>, 1 , 5001> dEmdz;
    Array<std::complex<double>, 1 , 5001> dEp2dz2;
    Array<std::complex<double>, 1 , 5001> dEm2dz2;
    Array<std::complex<double>, 1 , 5001> dEpdt0;
    Array<std::complex<double>, 1 , 5001> dEmdt0;
    Array<std::complex<double>, 1 , 5001> dEpdt;
    Array<std::complex<double>, 1 , 5001> dEmdt;
    Array<double, 1 , 5001> spectrum;
    double s = 501.;
    double k = 2;
    double max = 0;
    for(int x = 0; x < 5001; x++){
        Ep[x]= cos(2*M_PI*k/5001.*x*x/5001.);
        Em[x]= 0;
        spectrum[x] = 0.;
    }
    std::vector<double> t,p;
    Eigen::ArrayXf x = Eigen::ArrayXf::LinSpaced(5001, 0, 1);
    for(int i = 1; i < 5'000'000 ; i++)
    {
        Pp = Ep.abs().pow(2);
        Pm = Em.abs().pow(2);

        tmpdz[0] = Em[1]*sqrt(R1);
        tmpdz.rightCols(5001) = Ep;
        dEpdz =(tmpdz.middleCols(1,5001) - tmpdz.middleCols(0,5001))/dz;

        tmpdz[5001] = Ep[4999]*sqrt(R2);
        tmpdz.leftCols(5001) = Em;
        dEmdz =(tmpdz.middleCols(1,5001) - tmpdz.middleCols(0,5001))/dz;

        tmpdz2[0] = Em[2]*sqrt(R1);
        tmpdz2[1] = Em[1]*sqrt(R1);
        tmpdz2.middleCols(2,5001) = Ep;
        dEp2dz2 = (tmpdz2.middleCols(2, 5001) - 2.*tmpdz2.middleCols(1, 5001) + tmpdz2.middleCols(0,5001))/(dz*dz);

        tmpdz2[5002] = Ep[4998]*sqrt(R2);
        tmpdz2[5001] = Ep[4997]*sqrt(R2);
        tmpdz2.middleCols(0,5001) = Em;
        dEm2dz2 = (tmpdz2.middleCols(2, 5001) - 2.*tmpdz2.middleCols(1, 5001) + tmpdz2.middleCols(0,5001))/(dz*dz);

        dEpdt0 = -c/n*dEpdz;
        dEmdt0 =  c/n*dEmdz;

        dEpdt = c/n*(-dEpdz + ii*kpp/2.*(c/n)*(c/n)*dEp2dz2 - aw/2*Ep + 
                         -ii*gammaK*(Pp+2*Pm)*Ep + 
                         g0/2.*(Ep -1/Psat*(Pp+2*Pm)*Ep - T2*dEpdt0 + gc*T2*T2*(c/n)*(c/n)*dEp2dz2 + 
                     1/Psat*((2*T1+3*T2)*dEmdt0.conjugate()*Em*Ep + (T1+5/2.*T2)*Em.conjugate()*dEmdt0*Ep)));

        dEmdt = c/n*(dEmdz + ii*kpp/2.*(c/n)*(c/n)*dEm2dz2 - aw/2*Em + 
                     -ii*gammaK*(Pm+2*Pp)*Em + 
                         g0/2.*( Em -1/Psat*(Pm+2*Pp)*Em - T2*dEmdt0 + gc*T2*T2*(c/n)*(c/n)*dEm2dz2 + 
                     1/Psat*((2*T1+3*T2)*dEpdt0.conjugate()*Ep*Em + (T1+5/2.*T2)*Ep.conjugate()*dEpdt0*Em)));
        Ep = Ep + dEpdt*dt;
        Em = Em + dEmdt*dt;
        if(i % 200 == 0)
        {
            Eigen::Array<double, 5001, 1> y = (Pp + Pm).real();
            t.push_back(i);
            p.push_back(y.mean());
            max = std::max(max,y.mean());
            ui->leftPlot->yAxis->setRange(0, max);
            ui->leftPlot->xAxis->setRange(0, i);
            auto a = QVector<double>(t.data(), t.data() + t.size());
            auto b = QVector<double>(p.data(), p.data() + p.size());
            ui->leftPlot->graph(0)->setData(a,b);
            ui->leftPlot->replot();
            std::cout << i << " " << y.mean() << std::endl;


            double N = 5001;
            auto F  = J;
            Array<std::complex<double>, 1 ,5001> y2;
            fftw_complex *in, *out;
            fftw_plan p;
            in = (fftw_complex*) y.data();
            out = (fftw_complex*) y2.data();
            p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p); /* repeat as needed */
            fftw_destroy_plan(p);
            spectrum =  spectrum + y2.real()/y2.real().maxCoeff();
            auto a1 = QVector<double>(x.data(), x.data() + x.size());
            auto b1 = QVector<double>(spectrum.data(), spectrum.data() + spectrum.size());
            ui->rightPlot->yAxis->setRange(0,100);
            ui->rightPlot->xAxis->setRange(-.010,0.02);
            ui->rightPlot->graph(0)->setData(a1,b1);
            ui->rightPlot->replot();
            std::cout << y1 << std::endl;
            
        }
        

    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";






}

void printParams() {
    std::cout << "J: " << sim_params.J << std::endl;
    std::cout << "kpp: " << sim_params.kpp << std::endl;
    std::cout << "gammaK: " << sim_params.gammaK << std::endl;
    std::cout << "R1: " << sim_params.R1 << std::endl;
    std::cout << "R2: " << sim_params.R2 << std::endl;
    std::cout << "lm0: " << sim_params.lm0 << std::endl;
    std::cout << "Lc: " << sim_params.Lc << std::endl;
    std::cout << "df: " << sim_params.df << std::endl;
    std::cout << "mu: " << sim_params.mu << std::endl;
    std::cout << "T1: " << sim_params.T1 << std::endl;
    std::cout << "T2: " << sim_params.T2 << std::endl;
    std::cout << "n: " << sim_params.n << std::endl;
    std::cout << "Lmod: " << sim_params.Lmod << std::endl;
    std::cout << "aw: " << sim_params.aw << std::endl;
    std::cout << "Gamma: " << sim_params.Gamma << std::endl;
    std::cout << "h: " << sim_params.h << std::endl;
    std::cout << "dtperTr: " << sim_params.dtperTr << std::endl;
    std::cout << "numTr: " << sim_params.numTr << std::endl;
    std::cout << "Crnt: " << sim_params.Crnt << std::endl;
    std::cout << "gc: " << sim_params.gc << std::endl;
    std::cout << "gcfac: " << sim_params.gcfac << std::endl;
    std::cout << "lineshape: " << sim_params.ls << std::endl;
    std::cout << "useN2N3: " << sim_params.useN2N3 << std::endl;
    std::cout << "phi: " << sim_params.initphi << std::endl;
    std::cout << "plotprogress: " << sim_params.plotprogress << std::endl;
    std::cout << "plottheory: " << sim_params.plottheory << std::endl;
    std::cout << "shifttomin: " << sim_params.shifttomin << std::endl;
}

void MainWindow::on_Simulate_clicked()
{
    sim_params.J = ui->J->text().toDouble();
    sim_params.kpp = ui->Kpp->text().toDouble();
    sim_params.gammaK = ui->gammaK->text().toDouble();
    sim_params.R1 = ui->R1->text().toDouble();
    sim_params.R2 = ui->R2->text().toDouble();
    sim_params.lm0 = ui->lm0->text().toDouble();
    sim_params.Lc = ui->Lc->text().toDouble();
    sim_params.df = ui->df->text().toDouble();
    sim_params.mu = ui->mu->text().toDouble();
    sim_params.T1 = ui->T1->text().toDouble();
    sim_params.T2 = ui->T2->text().toDouble();
    sim_params.n = ui->n->text().toDouble();
    sim_params.Lmod = ui->Lmod->text().toDouble();
    sim_params.aw = ui->aw->text().toDouble();
    sim_params.Gamma = ui->Gamma->text().toDouble();
    sim_params.h = ui->h->text().toDouble();
    sim_params.dtperTr = ui->dtperTr->text().toInt();
    sim_params.numTr = ui->numTr->text().toInt();
    sim_params.Crnt = ui->Crnt->text().toDouble();
    sim_params.gc = ui->gc->isChecked();
    sim_params.gcfac = ui->gcfac->text().toDouble();
    sim_params.ls = static_cast<lineshape>(ui->Ls->currentIndex());
    sim_params.useN2N3 = ui->useN2N3->isChecked();
    sim_params.initphi = static_cast<phi>(ui->initphi->currentIndex());
    sim_params.plotprogress = ui->plotprogress->isChecked();
    sim_params.plottheory = ui->plottheory->isChecked();
    sim_params.shifttomin = ui->shifttomin->isChecked();
    printParams();
    simulate();
}
