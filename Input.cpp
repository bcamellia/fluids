#include "Input.h"

namespace InputManager
{
	namespace{
		short input_type;
		wiimote* remote;
	}
	void init(){
		input_type = INPUT_TYPE_KEYBOARD;
		remote = 0;
	}


	/*vmml::vec2i getMousePos(){
		return vmml::vec2i(sf::Mouse::getPosition().x,sf::Mouse::getPosition().y);
	}

	bool isMouseDownL(){
		return sf::Mouse::isButtonPressed(sf::Mouse::Button::Left);
	}

	bool isMouseDownR(){
		return sf::Mouse::isButtonPressed(sf::Mouse::Button::Right);
	}*/

	void on_state_change (wiimote &remote,
		state_change_flags  changed,
		const wiimote_state &new_state)
	{
		// we use this callback to set report types etc. to respond to key events
		//  (like the wiimote connecting or extensions (dis)connecting).

		// NOTE: don't access the public state from the 'remote' object here, as it will
		//		  be out-of-date (it's only updated via RefreshState() calls, and these
		//		  are reserved for the main application so it can be sure the values
		//		  stay consistent between calls).  Instead query 'new_state' only.

		// the wiimote just connected
		if(changed & CONNECTED)
		{
			// ask the wiimote to report everything (using the 'non-continous updates'
			//  default mode - updates will be frequent anyway due to the acceleration/IR
			//  values changing):

			// note1: you don't need to set a report type for Balance Boards - the
			//		   library does it automatically.

			// note2: for wiimotes, the report mode that includes the extension data
			//		   unfortunately only reports the 'BASIC' IR info (ie. no dot sizes),
			//		   so let's choose the best mode based on the extension status:
			if(new_state.ExtensionType != wiimote::BALANCE_BOARD)
			{
				if(new_state.bExtension)
					remote.SetReportType(wiimote::IN_BUTTONS_ACCEL_IR_EXT); // no IR dots
				else
					remote.SetReportType(wiimote::IN_BUTTONS_ACCEL_IR);		//    IR dots
			}
		}
		// a MotionPlus was detected
		if(changed & MOTIONPLUS_DETECTED)
		{
			// enable it if there isn't a normal extension plugged into it
			// (MotionPlus devices don't report like normal extensions until
			//  enabled - and then, other extensions attached to it will no longer be
			//  reported (so disable the M+ when you want to access them again).
			if(remote.ExtensionType == wiimote_state::NONE) {
				bool res = remote.EnableMotionPlus();
				_ASSERT(res);
			}
		}
		// an extension is connected to the MotionPlus
		else if(changed & MOTIONPLUS_EXTENSION_CONNECTED)
		{
			// We can't read it if the MotionPlus is currently enabled, so disable it:
			if(remote.MotionPlusEnabled())
				remote.DisableMotionPlus();
		}
		// an extension disconnected from the MotionPlus
		else if(changed & MOTIONPLUS_EXTENSION_DISCONNECTED)
		{
			// enable the MotionPlus data again:
			if(remote.MotionPlusConnected())
				remote.EnableMotionPlus();
		}
		// another extension was just connected:
		else if(changed & EXTENSION_CONNECTED)
		{
			// switch to a report mode that includes the extension data (we will
			//  loose the IR dot sizes)
			// note: there is no need to set report types for a Balance Board.
			if(!remote.IsBalanceBoard())
				remote.SetReportType(wiimote::IN_BUTTONS_ACCEL_IR_EXT);
		}
		// extension was just disconnected:
		else if(changed & EXTENSION_DISCONNECTED)
		{
			// use a non-extension report mode (this gives us back the IR dot sizes)
			remote.SetReportType(wiimote::IN_BUTTONS_ACCEL_IR);
		}
	}

	void setInputType(short type){
		input_type = type;
		switch(type){
		case INPUT_TYPE_WIIMOTE:
			findWiiRemote();
		}
	}

	void findWiiRemote(){
		remote = new wiimote;
		remote->ChangedCallback		= on_state_change;
		remote->CallbackTriggerFlags = (state_change_flags)(CONNECTED |
			EXTENSION_CHANGED |
			MOTIONPLUS_CHANGED);
		while(!remote->Connect(wiimote::FIRST_AVAILABLE)) {
		}
		remote->SetLEDs(IB8(1010));

	}

	void processMessages(){}
	/*	sf::Event msg;
		while(window.pollEvent(msg)){
			switch(msg.type){
			case sf::Event::MouseMoved:
				break;
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::KeyPressed:
				if(msg.key.code == sf::Keyboard::Escape){
					window.close();
				}
				break;
			case sf::Event::Resized:
				glViewport(0, 0, msg.size.width, msg.size.height);
				break;
			}
		}*/


	/*if(input_type == INPUT_TYPE_WIIMOTE){
	while(remote->RefreshState() == NO_CHANGE){
	Sleep(1); // // don't hog the CPU if nothing changed
	}
	}*/


	void destroy(){
		if(remote){
			remote->Disconnect();
		}
	}
}