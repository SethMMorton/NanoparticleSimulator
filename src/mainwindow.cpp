#include "mainwindow.h"
#include <QtDebug>
#include <QtGui>
#include <QtCore>
#include <QDir>
#include <QList>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    this->startPython();
    this->createWidgets();
    this->initUI();
    this->connectWidgets();

    // Set the default function to be sin
    this->functionField->setText("a*sin(b*x+c)+d");
    this->functionField->plotFunction();
    this->vbox[0]->changeParameter(0);
    this->vbox[1]->changeParameter(1);
    // Reset these values to 10 to correct a bug
    this->function->setMaxValue(0, 10.0);
    this->function->setMaxValue(1, 10.0);
}

void MainWindow::createWidgets()
{
    /* Creates all the widgets, i.e. run their constructors */

    // Load the python library module
    PyObject *pylib = this->loadPythonModule((char*) "pylib");

    // Now create the widgets
    this->plot          = new Plot(this);
    this->function      = new Function(pylib, this);
    this->functionField = new InputFunction(this);
    this->plotButton    = new QPushButton(tr("&Plot"), this);
    this->vbox[0]       = new ValueBox(tr("Parameter Value Picker 1"),
                                       0, this->function, this);
    this->vbox[1]       = new ValueBox(tr("Parameter Value Picker 2"),
                                       1, this->function, this);

    // Remove the reference count to this module
    Py_XDECREF(pylib);
}

void MainWindow::initUI()
{
    /* Sets up the layout of the window */

    // Define a central widget and a layout for the window
    this->setCentralWidget(new QWidget());
    this->mainLayout = new QVBoxLayout();
    this->setWindowTitle(tr("VisualFit"));

    // Set up the layout
    this->mainLayout->addWidget(this->vbox[0]);
    this->mainLayout->addWidget(this->vbox[1]);
    this->funcLay = new QHBoxLayout();
    this->funcLay->addWidget(new QLabel(tr("Enter a function: f(x) = ")));
    this->funcLay->addWidget(this->functionField);
    this->funcLay->addWidget(this->plotButton);
    this->mainLayout->addLayout(this->funcLay);
    this->mainLayout->addWidget(this->plot);

    // Add the widgets to the central widget
    this->centralWidget()->setLayout(this->mainLayout);
}

void MainWindow::connectWidgets()
{
    /* Connect the main widgets in this window */

    // If the plot button is pressed (or return is pressed), plot the function
    connect(this->plotButton,    SIGNAL(clicked()),
            this->functionField, SLOT(plotFunction()));
    connect(this->functionField, SIGNAL(returnPressed()),
            this->functionField, SLOT(plotFunction()));

    // If the function is broadcasted, have the function class catch and  parse it
    connect(this->functionField, SIGNAL(broadcastFunction(QString)),
            this->function,      SLOT(parseFunction(QString)));

    // When the function class emits the function, have the plot grab it
    connect(this->function, SIGNAL(hasPlotableFunction(PyObject*)),
            this->plot,     SLOT(plotFunction(PyObject*)));
}

MainWindow::~MainWindow()
{

}
