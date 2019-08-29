#include "simple_time.h"

namespace simple {
	namespace time {
		
		void clock::SetFps(float _fps) {
			fps = _fps;
			time_step = 1.0f / fps;
		}

		float clock::Update() {
			accumulator += GetDeltaTime();
			float out = accumulator;
			if (accumulator >= time_step || fps == -1)
			{
				accumulator = 0.0f; return out;
			}
			else
				return 0.0f;
		}

		float clock::GetDeltaTime() {
			time_current = sys_clock::now();
			std::chrono::duration<float> elapsed_time = time_current - time_previous;
			time_previous = time_current;
			return elapsed_time.count();
		}
	}
}