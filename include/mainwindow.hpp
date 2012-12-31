#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui/QMainWindow>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QWidget>
#include <QtGui/QPushButton>
//#include "plot.h"
//#include "valuebox.h"
//#include "function.h"
//#include "inputfunction.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    void createWidgets();
    void initUI();
    void connectWidgets();
    void startPython();
    PyObject *loadPythonModule(const char* modname);

private:
    QVBoxLayout *mainLayout;
    QHBoxLayout *funcLay;
    Plot *plot;
    ValueBox *vbox[2];
    Function *function;
    InputFunction *functionField;
    QPushButton *plotButton;

};

#endif // MAINWINDOW_H
