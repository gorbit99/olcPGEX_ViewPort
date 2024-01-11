#pragma once

#include "olcPixelGameEngine.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

// Declarations
namespace olc {
    class ViewPort : public olc::PGEX {
    public:
        ViewPort();
        ViewPort(std::vector<vf2d> vertices, vf2d offset = {0, 0});
        virtual ~ViewPort();
        void addPoint(vf2d point);
        void clear();
        void drawEdges();
        void setOffset(vf2d offset);

        static ViewPort rectViewPort(vf2d topLeft,
                                     vf2d size,
                                     olc::vf2d offset = {0, 0});

        void DrawDecal(const olc::vf2d &pos,
                       olc::Decal *decal,
                       const olc::vf2d &scale = {1.0f, 1.0f},
                       const olc::Pixel &tint = olc::WHITE) const;
        void DrawPartialDecal(const olc::vf2d &pos,
                              olc::Decal *decal,
                              const olc::vf2d &source_pos,
                              const olc::vf2d &source_size,
                              const olc::vf2d &scale = {1.0f, 1.0f},
                              const olc::Pixel &tint = olc::WHITE) const;
        void DrawPartialDecal(const vf2d &pos,
                              const vf2d &size,
                              Decal *decal,
                              const vf2d source_pos,
                              const vf2d &source_size,
                              const Pixel &tint = olc::WHITE) const;
        void DrawExplicitDecal(olc::Decal *decal,
                               const olc::vf2d *pos,
                               const olc::vf2d *uv,
                               const olc::Pixel *col,
                               uint32_t elements = 4) const;
        void DrawWarpedDecal(Decal *decal,
                             const vf2d (&pos)[4],
                             const Pixel &tint = WHITE) const;
        void DrawWarpedDecal(Decal *decal,
                             const vf2d *pos,
                             const Pixel &tint = WHITE) const;
        void DrawWarpedDecal(Decal *decal,
                             const std::array<vf2d, 4> &pos,
                             const Pixel &tint = WHITE) const;
        void DrawPartialWarpedDecal(Decal *decal,
                                    const vf2d (&pos)[4],
                                    const vf2d &source_pos,
                                    const vf2d &source_size,
                                    const Pixel &tint = WHITE) const;
        void DrawPartialWarpedDecal(Decal *decal,
                                    const vf2d *pos,
                                    const vf2d &source_pos,
                                    const vf2d &source_size,
                                    const Pixel &tint = WHITE) const;
        void DrawPartialWarpedDecal(Decal *decal,
                                    const std::array<vf2d, 4> &pos,
                                    const vf2d &source_pos,
                                    const vf2d &source_size,
                                    const Pixel &tint = WHITE) const;
        void DrawRotatedDecal(const vf2d &pos,
                              Decal *decal,
                              const float fAngle,
                              const vf2d &center = {0.0f, 0.0f},
                              const vf2d &scale = {1.0f, 1.0f},
                              const Pixel &tint = WHITE) const;
        void DrawPartialRotatedDecal(const vf2d &pos,
                                     Decal *decal,
                                     const float fAngle,
                                     const vf2d &center,
                                     const vf2d &source_pos,
                                     const vf2d &source_size,
                                     const vf2d &scale = {1.0f, 1.0f},
                                     const Pixel &tint = WHITE) const;
        void DrawRectDecal(const vf2d &pos,
                           const vf2d &size,
                           const Pixel col = WHITE) const;
        void FillRectDecal(const vf2d &pos,
                           const vf2d &size,
                           const Pixel col = WHITE) const;
        void GradientFillRectDecal(const vf2d &pos,
                                   const vf2d &size,
                                   const Pixel colTL,
                                   const Pixel colBL,
                                   const Pixel colBR,
                                   const Pixel colTR) const;
        void DrawPolygonDecal(Decal *decal,
                              const std::vector<vf2d> &pos,
                              const std::vector<vf2d> &uv,
                              const Pixel tint = WHITE) const;
        void DrawPolygonDecal(Decal *decal,
                              const std::vector<vf2d> &pos,
                              const std::vector<float> &depth,
                              const std::vector<vf2d> &uv,
                              const Pixel tint = WHITE) const;
        void DrawPolygonDecal(Decal *decal,
                              const std::vector<vf2d> &pos,
                              const std::vector<vf2d> &uv,
                              const std::vector<Pixel> &tint) const;
        void DrawLineDecal(const vf2d &pos1,
                           const vf2d &pos2,
                           Pixel p = WHITE) const;

    private:
        void drawClippedDecal(Decal *decal,
                              const vf2d *points,
                              const vf2d *uvs,
                              const Pixel *col,
                              uint32_t elements = 0) const;
        void drawClippedPolygonDecal(Decal *decal,
                              const vf2d *points,
                              const vf2d *uvs,
                              const float *depth,
                              const Pixel tint,
                              uint32_t elements = 0) const;
        
        static bool ccw(vf2d A,vf2d B,vf2d C);
        static bool intersect(vf2d A,vf2d B,vf2d C,vf2d D);
        static float lineSegmentIntersect(vf2d lineA,
                                          vf2d lineB,
                                          vf2d segmentA,
                                          vf2d segmentB);
        static float directionFromLine(vf2d lineA, vf2d lineB, vf2d point);

        std::vector<vf2d> clipVertices;
        olc::vf2d offset;
    };
} // namespace olc

// Definitions

#ifdef OLC_PGEX_VIEWPORT
#undef OLC_PGEX_VIEWPORT

olc::ViewPort::ViewPort() {
}
olc::ViewPort::~ViewPort() {
}

olc::ViewPort::ViewPort(std::vector<vf2d> vertices, olc::vf2d offset)
        : clipVertices{vertices},
          offset{offset} {
}

void olc::ViewPort::addPoint(vf2d point) {
    clipVertices.push_back(point);
}

void olc::ViewPort::clear() {
    clipVertices.clear();
}

void olc::ViewPort::drawEdges() {
    for (auto i = 0u; i < clipVertices.size(); i++) {
        auto current = clipVertices[i] + offset;
        auto next = clipVertices[(i + 1) % clipVertices.size()] + offset;

        pge->DrawLineDecal(current, next, olc::RED);
    }
}

void olc::ViewPort::setOffset(vf2d offset) {
    this->offset = offset;
}

olc::ViewPort
        olc::ViewPort::rectViewPort(vf2d topLeft, vf2d size, olc::vf2d offset) {
    return {{
                    topLeft,
                    {topLeft.x, topLeft.y + size.y},
                    topLeft + size,
                    {topLeft.x + size.x, topLeft.y},
            },
            offset};
}

void olc::ViewPort::DrawDecal(const olc::vf2d &pos,
                              olc::Decal *decal,
                              const olc::vf2d &scale,
                              const olc::Pixel &tint) const {
    std::vector<olc::vf2d> points{
            pos,
            {pos.x, pos.y + decal->sprite->height * scale.y},
            {pos.x + decal->sprite->width * scale.x,
             pos.y + decal->sprite->height * scale.y},
            {pos.x + decal->sprite->width * scale.x, pos.y},
    };
    DrawWarpedDecal(decal, points.data(), tint);
}

void olc::ViewPort::DrawPartialDecal(const olc::vf2d &pos,
                                     olc::Decal *decal,
                                     const olc::vf2d &source_pos,
                                     const olc::vf2d &source_size,
                                     const olc::vf2d &scale,
                                     const olc::Pixel &tint) const {
    DrawPartialDecal(pos, source_size * scale, decal, source_pos, source_size, tint);
}

void olc::ViewPort::DrawPartialDecal(const vf2d &pos,
                                     const vf2d &size,
                                     Decal *decal,
                                     const vf2d source_pos,
                                     const vf2d &source_size,
                                     const Pixel &tint) const {
    std::vector<vf2d> points{
            pos,
            {pos.x, pos.y + size.y},
            pos + size,
            {pos.x + size.x, pos.y},
    };
    DrawPartialWarpedDecal(decal, points.data(), source_pos, source_size, tint);
}

void olc::ViewPort::DrawExplicitDecal(olc::Decal *decal,
                                      const olc::vf2d *pos,
                                      const olc::vf2d *uv,
                                      const olc::Pixel *col,
                                      uint32_t elements) const {
    drawClippedDecal(decal, pos, uv, col, elements);
}

void olc::ViewPort::DrawWarpedDecal(Decal *decal,
                                    const vf2d (&pos)[4],
                                    const Pixel &tint) const {
    DrawWarpedDecal(decal, (const vf2d *)pos, tint);
}
void olc::ViewPort::DrawWarpedDecal(Decal *decal,
                                    const vf2d *pos,
                                    const Pixel &tint) const {
    std::vector<float> w{ 1, 1, 1, 1 };
    std::vector<olc::vf2d> newPos;
	newPos.resize(4);
    std::vector<vf2d> uvs{
            {0, 0},
            {0, 1},
            {1, 1},
            {1, 0},
    };
    std::vector<Pixel> cols{
            tint,
            tint,
            tint,
            tint,
    };

    olc::vf2d vInvScreenSize={ 1.0f / pge->GetScreenSize().x, 1.0f / pge->GetScreenSize().y };

    olc::vf2d center;
	float rd = ((pos[2].x - pos[0].x) * (pos[3].y - pos[1].y) - (pos[3].x - pos[1].x) * (pos[2].y - pos[0].y));
	if (rd != 0)
	{
		rd = 1.0f / rd;
		float rn = ((pos[3].x - pos[1].x) * (pos[0].y - pos[1].y) - (pos[3].y - pos[1].y) * (pos[0].x - pos[1].x)) * rd;
		float sn = ((pos[2].x - pos[0].x) * (pos[0].y - pos[1].y) - (pos[2].y - pos[0].y) * (pos[0].x - pos[1].x)) * rd;
		if (!(rn < 0.f || rn > 1.f || sn < 0.f || sn > 1.f)) center = pos[0] + rn * (pos[2] - pos[0]);
		float d[4];	for (int i = 0; i < 4; i++)	d[i] = (pos[i] - center).mag();
		for (int i = 0; i < 4; i++)
		{
			float q = d[i] == 0.0f ? 1.0f : (d[i] + d[(i + 2) & 3]) / d[(i + 2) & 3];
			uvs[i] *= q; w[i] *= q;
		}

        drawClippedPolygonDecal(decal, pos, uvs.data(), w.data(), tint, 4);
    }
}
void olc::ViewPort::DrawWarpedDecal(Decal *decal,
                                    const std::array<vf2d, 4> &pos,
                                    const Pixel &tint) const {
    DrawWarpedDecal(decal, pos.data(), tint);
}

void olc::ViewPort::DrawPartialWarpedDecal(Decal *decal,
                                           const vf2d (&pos)[4],
                                           const vf2d &source_pos,
                                           const vf2d &source_size,
                                           const Pixel &tint) const {
    DrawPartialWarpedDecal(decal,
                           (const vf2d *)pos,
                           source_pos,
                           source_size,
                           tint);
}

void olc::ViewPort::DrawPartialWarpedDecal(Decal *decal,
                                           const vf2d *pos,
                                           const vf2d &source_pos,
                                           const vf2d &source_size,
                                           const Pixel &tint) const {
    olc::vf2d sourceUvPos =
            source_pos
            / olc::vf2d{static_cast<float>(decal->sprite->width),
                        static_cast<float>(decal->sprite->height)};
    olc::vf2d sourceUvSize =
            source_size
            / olc::vf2d{static_cast<float>(decal->sprite->width),
                        static_cast<float>(decal->sprite->height)};
    std::vector<vf2d> uvs{
            sourceUvPos,
            {sourceUvPos.x, sourceUvPos.y + sourceUvSize.y},
            sourceUvPos + sourceUvSize,
            {sourceUvPos.x + sourceUvSize.x, sourceUvPos.y},
    };
    std::vector<Pixel> cols{
            tint,
            tint,
            tint,
            tint,
    };

    std::vector<float>ws{1,1,1,1};

    olc::vf2d center;
	float rd = ((pos[2].x - pos[0].x) * (pos[3].y - pos[1].y) - (pos[3].x - pos[1].x) * (pos[2].y - pos[0].y));
	if (rd != 0)
	{
		rd = 1.0f / rd;
		float rn = ((pos[3].x - pos[1].x) * (pos[0].y - pos[1].y) - (pos[3].y - pos[1].y) * (pos[0].x - pos[1].x)) * rd;
		float sn = ((pos[2].x - pos[0].x) * (pos[0].y - pos[1].y) - (pos[2].y - pos[0].y) * (pos[0].x - pos[1].x)) * rd;
		if (!(rn < 0.f || rn > 1.f || sn < 0.f || sn > 1.f)) center = pos[0] + rn * (pos[2] - pos[0]);
		float d[4];	for (int i = 0; i < 4; i++)	d[i] = (pos[i] - center).mag();
		for (int i = 0; i < 4; i++)
		{
			float q = d[i] == 0.0f ? 1.0f : (d[i] + d[(i + 2) & 3]) / d[(i + 2) & 3];
			uvs[i] *= q; ws[i] *= q;
		}

        drawClippedPolygonDecal(decal, pos, uvs.data(), ws.data(), tint, 4);
    }
}

void olc::ViewPort::DrawPartialWarpedDecal(Decal *decal,
                                           const std::array<vf2d, 4> &pos,
                                           const vf2d &source_pos,
                                           const vf2d &source_size,
                                           const Pixel &tint) const {
    DrawPartialWarpedDecal(decal, pos.data(), source_pos, source_size, tint);
}

void olc::ViewPort::DrawRotatedDecal(const vf2d &pos,
                                     Decal *decal,
                                     const float fAngle,
                                     const vf2d &center,
                                     const vf2d &scale,
                                     const Pixel &tint) const {
    auto sin = std::sin(fAngle);
    auto cos = std::cos(fAngle);

    std::vector<vf2d> points{
            -center * scale,
            olc::vf2d{-center.x, decal->sprite->height - center.y} * scale,
            olc::vf2d{decal->sprite->width - center.x,
                      decal->sprite->height - center.y}
                    * scale,
            olc::vf2d{decal->sprite->width - center.x, -center.y} * scale,
    };

    for (auto i = 0u; i < points.size(); i++) {
        points[i] = pos
                    + olc::vf2d{points[i].x * cos - points[i].y * sin,
                                points[i].x * sin + points[i].y * cos};
    }

    DrawWarpedDecal(decal, points.data(), tint);
}

void olc::ViewPort::DrawPartialRotatedDecal(const vf2d &pos,
                                            Decal *decal,
                                            const float fAngle,
                                            const vf2d &center,
                                            const vf2d &source_pos,
                                            const vf2d &source_size,
                                            const vf2d &scale,
                                            const Pixel &tint) const {
    auto sin = std::sin(fAngle);
    auto cos = std::cos(fAngle);

    std::vector<vf2d> points{
            -center * scale,
            olc::vf2d{-center.x, source_size.y - center.y} * scale,
            (source_size - center) * scale,
            olc::vf2d{source_size.x - center.x, -center.y} * scale,
    };

    for (auto i = 0u; i < points.size(); i++) {
        points[i] = pos
                    + olc::vf2d{points[i].x * cos - points[i].y * sin,
                                points[i].x * sin + points[i].y * cos};
    }

    DrawPartialWarpedDecal(decal, points.data(), source_pos, source_size, tint);
}

void olc::ViewPort::DrawRectDecal(const vf2d &pos,
                                  const vf2d &size,
                                  const Pixel col) const {
    std::vector<vf2d> points{
            pos,
            {pos.x, pos.y + size.y},
            pos + size,
            {pos.x + size.x, pos.y},
    };

    // Ideally we use the wireframe mode just like the PGE, 
    // however we can't save the current decal mode which 
    // can impact some applications so instead we draw 4
    // lines.

    DrawLineDecal(points[0],points[1],col);
    DrawLineDecal(points[1],points[2],col);
    DrawLineDecal(points[2],points[3],col);
    DrawLineDecal(points[3],points[0],col);
}

void olc::ViewPort::FillRectDecal(const vf2d &pos,
                                  const vf2d &size,
                                  const Pixel col) const {
    std::vector<vf2d> points{
            pos,
            {pos.x, pos.y + size.y},
            pos + size,
            {pos.x + size.x, pos.y},
    };
    std::vector<vf2d> uvs{
            {0, 0},
            {0, 1},
            {1, 1},
            {1, 0},
    };

    DrawPolygonDecal(nullptr, points, uvs, col);
}

void olc::ViewPort::GradientFillRectDecal(const vf2d &pos,
                                          const vf2d &size,
                                          const Pixel colTL,
                                          const Pixel colBL,
                                          const Pixel colBR,
                                          const Pixel colTR) const {
    std::vector<vf2d> points{
            pos,
            {pos.x, pos.y + size.y},
            pos + size,
            {pos.x + size.x, pos.y},
    };

    std::vector<vf2d> uvs{
            {0, 0},
            {0, 1},
            {1, 1},
            {1, 0},
    };

    std::vector<Pixel> colors{
            colTL,
            colBL,
            colBR,
            colTR,
    };

    drawClippedDecal(nullptr, points.data(), uvs.data(), colors.data(), points.size());
}

void olc::ViewPort::DrawPolygonDecal(Decal *decal,
                                     const std::vector<vf2d> &pos,
                                     const std::vector<vf2d> &uv,
                                     const Pixel tint) const {
    std::vector<Pixel> colors;
    colors.resize(pos.size());
    for (auto i = 0u; i < colors.size(); i++) {
        colors[i] = tint;
    }

    drawClippedDecal(decal, pos.data(), uv.data(), colors.data(), pos.size());
}

void olc::ViewPort::DrawPolygonDecal(Decal *decal,
                                     const std::vector<vf2d> &pos,
                                     const std::vector<float> &depth,
                                     const std::vector<vf2d> &uv,
                                     const Pixel tint) const {
   drawClippedPolygonDecal(decal, pos.data(), uv.data(), depth.data(), tint, pos.size());
}

void olc::ViewPort::DrawPolygonDecal(Decal *decal,
                                     const std::vector<vf2d> &pos,
                                     const std::vector<vf2d> &uv,
                                     const std::vector<Pixel> &tint) const {
    drawClippedDecal(decal, pos.data(), uv.data(), tint.data(), pos.size());
}

void olc::ViewPort::DrawLineDecal(const vf2d &pos1,
                                  const vf2d &pos2,
                                  Pixel p) const {
    vf2d posA = pos1 + offset;
    vf2d posB = pos2 + offset;

    for (auto i = 0u; i < clipVertices.size(); i++) {
        auto clipA = clipVertices[i] + offset;
        auto clipB = clipVertices[(i + 1) % clipVertices.size()] + offset;

        auto intersection = lineSegmentIntersect(clipA, clipB, posA, posB);
        if (intersection < 0 || intersection > 1) {
            continue;
        }

        auto clipDirection = directionFromLine(clipA, clipB, posA);
        auto intersectionPoint = posA + (posB - posA) * intersection;

        if (clipDirection >= 0) {
            posA = intersectionPoint;
        } else {
            posB = intersectionPoint;
        }
    }
    

    // Inside check. Draw a ray to the edge of the screen and count the times
    // it intersects. When odd, we are inside a shape, when even we are outside
    // of it.

    vf2d leftEdgeA = {0.f,posA.y};
    vf2d leftEdgeB = {0.f,posB.y};

    int leftEdgeIntersectionsA = 0;
    int leftEdgeIntersectionsB = 0;
    for (auto i = 0u; i < clipVertices.size(); i++) {
        auto clipA = clipVertices[i] + offset;
        auto clipB = clipVertices[(i + 1) % clipVertices.size()] + offset;
        auto leftEdgeIntersectA = intersect(clipA, clipB, leftEdgeA, posA);
        auto leftEdgeIntersectB = intersect(clipA, clipB, leftEdgeB, posB);

        if (leftEdgeIntersectA) {
           leftEdgeIntersectionsA++;
        }
        if (leftEdgeIntersectB) {
           leftEdgeIntersectionsB++;
        }
    }

    // If we found an intersection, we are drawing this line.
    // 
    // Otherwise, if either count is odd, one point is at 
    // least inside the shape, so render it.
    if (leftEdgeIntersectionsA % 2 == 1 || leftEdgeIntersectionsB % 2 == 1) {
        pge->DrawLineDecal(posA, posB, p);
    }
}

void olc::ViewPort::drawClippedDecal(Decal *decal,
                                     const vf2d *points,
                                     const vf2d *uvs,
                                     const Pixel *col,
                                     uint32_t elements) const {
    std::vector<vf2d> outputList{points, points + elements};
    std::vector<vf2d> outputUvs{uvs, uvs + elements};
    std::vector<Pixel> outputCols{col, col + elements};

    for (auto i = 0u; i < clipVertices.size(); i++) {
        auto clipA = clipVertices[i];
        auto clipB = clipVertices[(i + 1) % 4];

        auto inputList{outputList};
        auto inputUvs{outputUvs};
        auto inputCols{outputCols};
        outputList.clear();
        outputUvs.clear();
        outputCols.clear();

        for (auto i = 0u; i < inputList.size(); i++) {
            auto polygonA = inputList[i];
            auto polygonB = inputList[(i + 1) % inputList.size()];
            auto uvA = inputUvs[i];
            auto uvB = inputUvs[(i + 1) % inputList.size()];
            auto colA = inputCols[i];
            auto colB = inputCols[(i + 1) % inputList.size()];

            auto intersection =
                    lineSegmentIntersect(clipA, clipB, polygonA, polygonB);
            auto intersectionPoint =
                    polygonA + (polygonB - polygonA) * intersection;
            auto intersectionUv = uvA + (uvB - uvA) * intersection;
            auto intersectionCol = PixelLerp(colA, colB, intersection);

            float aDirection = directionFromLine(clipA, clipB, polygonA);
            float bDirection = directionFromLine(clipA, clipB, polygonB);

            if (bDirection <= 0) {
                if (aDirection > 0) {
                    outputList.push_back(intersectionPoint);
                    outputUvs.push_back(intersectionUv);
                    outputCols.push_back(intersectionCol);
                }
                outputList.push_back(polygonB);
                outputUvs.push_back(uvB);
                outputCols.push_back(colB);
            } else if (aDirection <= 0) {
                outputList.push_back(intersectionPoint);
                outputUvs.push_back(intersectionUv);
                outputCols.push_back(intersectionCol);
            }
        }
    }

    if (outputList.size() == 0) {
        return;
    }

    for (auto &point : outputList) {
        point += offset;
    }

    pge->DrawExplicitDecal(decal,
                           outputList.data(),
                           outputUvs.data(),
                           outputCols.data(),
                           outputList.size());
}
void olc::ViewPort::drawClippedPolygonDecal(Decal *decal,
                                     const vf2d *points,
                                     const vf2d *uvs,
                                     const float *depth,
                                     const Pixel tint,
                                     uint32_t elements) const {
    std::vector<vf2d> outputList{points, points + elements};
    std::vector<vf2d> outputUvs{uvs, uvs + elements};
    std::vector<float> outputDepths{depth, depth + elements};

    for (auto i = 0u; i < clipVertices.size(); i++) {
        auto clipA = clipVertices[i];
        auto clipB = clipVertices[(i + 1) % 4];

        auto inputList{outputList};
        auto inputUvs{outputUvs};
        auto inputWs{outputDepths};
        outputList.clear();
        outputUvs.clear();
        outputDepths.clear();

        for (auto i = 0u; i < inputList.size(); i++) {
            auto polygonA = inputList[i];
            auto polygonB = inputList[(i + 1) % inputList.size()];
            auto uvA = inputUvs[i];
            auto uvB = inputUvs[(i + 1) % inputList.size()];
            auto Wa = inputWs[i];
            auto Wb = inputWs[(i + 1) % inputList.size()];

            auto intersection =
                    lineSegmentIntersect(clipA, clipB, polygonA, polygonB);
            auto intersectionPoint =
                    polygonA + (polygonB - polygonA) * intersection;
            auto intersectionUv = uvA + (uvB - uvA) * intersection;
            auto intersectionDepth = Wa + (Wb - Wa) * intersection;

            float aDirection = directionFromLine(clipA, clipB, polygonA);
            float bDirection = directionFromLine(clipA, clipB, polygonB);

            if (bDirection <= 0) {
                if (aDirection > 0) {
                    outputList.push_back(intersectionPoint);
                    outputUvs.push_back(intersectionUv);
                    outputDepths.push_back(intersectionDepth);
                }
                outputList.push_back(polygonB);
                outputUvs.push_back(uvB);
                outputDepths.push_back(Wb);
            } else if (aDirection <= 0) {
                outputList.push_back(intersectionPoint);
                outputUvs.push_back(intersectionUv);
                outputDepths.push_back(intersectionDepth);
            }
        }
    }

    for (auto &point : outputList) {
        point += offset;
    }

    pge->DrawPolygonDecal(decal,
                          outputList,
                          outputDepths,
                          outputUvs,
                          tint);
}

bool olc::ViewPort::ccw(vf2d A,vf2d B,vf2d C) {
    return (C.y-A.y) * (B.x-A.x) > (B.y-A.y) * (C.x-A.x);
}

bool olc::ViewPort::intersect(vf2d A,vf2d B,vf2d C,vf2d D) {
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D);
}

float olc::ViewPort::lineSegmentIntersect(vf2d lineA,
                                          vf2d lineB,
                                          vf2d segmentA,
                                          vf2d segmentB) {
    return ((lineA.x - segmentA.x) * (lineA.y - lineB.y)
            - (lineA.y - segmentA.y) * (lineA.x - lineB.x))
           / ((lineA.x - lineB.x) * (segmentA.y - segmentB.y)
              - (lineA.y - lineB.y) * (segmentA.x - segmentB.x));
}

float olc::ViewPort::directionFromLine(vf2d lineA, vf2d lineB, vf2d point) {
    return (lineB.x - lineA.x) * (point.y - lineA.y)
           - (point.x - lineA.x) * (lineB.y - lineA.y);
}

#endif