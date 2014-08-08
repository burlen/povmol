#include <QApplication>
#include "PovmolMainWindow.h"

int main( int argc, char** argv )
{
  QApplication app(argc, argv);

  PovmolMainWindow win;
  win.show();

  return app.exec();
}
