#pragma once
#include "simple_math.h"

namespace simple {

	namespace gui {
		constexpr int max_widget_count = 256;

		struct widget {
			char* id;
			float height;

			bool is_hot;
			bool is_active;
		};

		struct string : public widget {
			string(char* _id) { id = _id; }
		};

		struct button : public widget {

		};

		struct menu_break : public widget {

		};

		struct menu_end : public widget {

		};

		class menu {
		private:
			math::vec2 position;
			math::vec2 mouse_position;

			widget* widgets[max_widget_count];
			int count;

			float padding = 16.0f;
			float spacing = 8.0f;
			float indent = 0.0f;

		public:
			void SetPadding(const float& f) { padding = f; }
			void SetSpacing(const float& f) { spacing = f; }
			void SetIndent(const float& f) { spacing = f; }

			void Start(const float& x, const float& y) {
				position = { x, y };
			}

			void GetMouseStatus(const float& mouse_x, const float& mouse_y) {
				mouse_position = { mouse_x, mouse_y };
			}

			void AddString(char* id, const char* value) {

				widgets[count] = new string(id);
			}

			void End() {
				for (auto w : widgets) {
					// Update hot and active
				}
				delete[] widgets;
			}
		};
	}
}