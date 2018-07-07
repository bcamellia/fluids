#ifndef INPUTCONTROL_H
#define INPUTCONTROL_H

#include <wiimote.h>
#include <vector.hpp>
#include "conversion.h"

#define INPUT_TYPE_KEYBOARD 0
#define INPUT_TYPE_WIIMOTE 1

namespace InputManager{
	void init();
	void setInputType(short type);
	/*vmml::vec2i getMousePos();
	bool isMouseDownL();
	bool isMouseDownR();*/
	void findWiiRemote();
	void destroy();
	void processMessages();
};

#endif