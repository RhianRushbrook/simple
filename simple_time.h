#pragma once
#include <chrono>

namespace simple {

	namespace time {

		typedef std::chrono::system_clock sys_clock;
		typedef std::chrono::time_point<std::chrono::system_clock> time_point;

		class clock {
		private:
			float fps = 60.0f;
			float time_step = 1.0f / fps;
			float accumulator = 0.0f;

			time_point time_previous = sys_clock::now();
			time_point time_current = sys_clock::now();

		public:
			void SetFps(float _fps);
			float GetDeltaTime();
			float Update();
		};
	}
}
