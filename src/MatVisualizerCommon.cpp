#include "MatVisualizerCommon.h"

#ifdef QMAT_WITH_POLYSCOPE

#include <cmath>

std::array<float,3> HsvToRgb(float h, float s, float v)
{
    float r = 0, g = 0, b = 0;
    int   i = (int)(h * 6.0f);
    float f = h * 6.0f - i;
    float p = v * (1.0f - s);
    float q = v * (1.0f - f * s);
    float t = v * (1.0f - (1.0f - f) * s);
    switch (i % 6) {
        case 0: r=v; g=t; b=p; break;
        case 1: r=q; g=v; b=p; break;
        case 2: r=p; g=v; b=t; break;
        case 3: r=p; g=q; b=v; break;
        case 4: r=t; g=p; b=v; break;
        case 5: r=v; g=p; b=q; break;
    }
    return {r, g, b};
}

std::array<float,3> StructIdColor(int struct_id)
{
    if (struct_id < 0) return {0.5f, 0.5f, 0.5f};
    const float golden = 0.618033988749895f;
    float hue = std::fmod(struct_id * golden, 1.0f);
    return HsvToRgb(hue, 0.85f, 0.95f);
}

#endif  // QMAT_WITH_POLYSCOPE
