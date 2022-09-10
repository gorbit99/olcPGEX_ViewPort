#include <vector>
#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#define OLC_PGEX_VIEWPORT
#include "olcPGEX_ViewPort.h"

class Example : public olc::PixelGameEngine {
public:
    olc::ViewPort viewPort;
    olc::Renderable uvChecker;
    olc::vf2d pos{10, 10};
    float angle = 0;

    bool OnUserCreate() override {
        sAppName = "OLC PGEX ViewPort";
        uvChecker.Load("testImage.png");
        viewPort = olc::ViewPort::rectViewPort({20, 20}, {100, 100});
        viewPort.setOffset({100, 100});
        return true;
    }

    bool OnUserUpdate(float fElapsedTime) override {
        if (GetKey(olc::A).bHeld) {
            pos.x -= fElapsedTime * 100;
        }
        if (GetKey(olc::D).bHeld) {
            pos.x += fElapsedTime * 100;
        }
        if (GetKey(olc::W).bHeld) {
            pos.y -= fElapsedTime * 100;
        }
        if (GetKey(olc::S).bHeld) {
            pos.y += fElapsedTime * 100;
        }
        if (GetKey(olc::Q).bHeld) {
            angle -= fElapsedTime * 3.141592 / 2;
        }
        if (GetKey(olc::E).bHeld) {
            angle += fElapsedTime * 3.141592 / 2;
        }

        viewPort.DrawRotatedDecal(pos,
                                  uvChecker.Decal(),
                                  angle,
                                  {uvChecker.Sprite()->width / 2.0f,
                                   uvChecker.Sprite()->height / 2.0f},
                                  {0.1f, 0.1f});

        viewPort.drawEdges();
        return true;
    }
};

int main() {
    Example example;
    if (example.Construct(640, 480, 1, 1)) {
        example.Start();
    }

    return 0;
}
