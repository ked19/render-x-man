#include "main.h"

void main()
{
	pMain_win = make_window();
	pMain_win->show();
	pRender_glWin->show();

	Fl::run();
}