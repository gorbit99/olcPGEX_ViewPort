# olcPGEX_ViewPort

An extension that allows you to draw decals only to a certain part of the
screen.

# Usage

Include the header file, and only once, somewhere before that include in a
source file, include the following line:

```cpp
#define OLC_PGEX_VIEWPORT
```

Afterwards, you can create a viewport object in one of the following ways:

```cpp
//Creates an empty viewport, and fills it up
olc::ViewPort viewport = olc::ViewPort();
addPoint({10, 10});
addPoint({10, 100});
addPoint({100, 100});

//Creates a viewport by specifying its points
//Make sure the points are in counter-clockwise order
std::vector<olc::vf2d> viewportPoints = {
        {20, 20},
        {20, 100},
        {100, 100},
        {100, 20},
    }
olc::ViewPort viewport = olc::ViewPort(viewportPoints);

//Creates a rectangular viewport
olc::ViewPort viewport = olc::ViewPort::rectViewPort(topLeft, size);
```

Afterwards, you can draw clipped decals, by calling the corresponding method on
the viewport object:

```cpp
viewport.DrawDecal(myDecal, position);
```

Besides text drawing, all decal methods are supported in the same way as in the
Pixel Game Engine

You can also set the offset of the viewport in one of the following ways:

```cpp
//Pass it into the constructor
olc::ViewPort viewport = olc::ViewPort(points, offset);

//Set using method
viewport.setOffset(offset);
```

Afterwards, both the drawings and the visible area will be shifted by that
amount. Useful if you want to add an offset to the drawings, like if it was a
window that's being moved.

To debug the viewport, call the `drawEdges()` method on it. This will outline
the shape so you know exactly where it is.

# Performance

Performance is dependent on how many vertices make up both the decal and the
viewport. With a rectangular viewport, and a decal with 4 vertices (most cases),
the performance impact is basically unnoticable.

# Limitations

Like in the Pixel Game Engine, drawing concave decals or using concave viewports 
is undefined behaviour. It might work, it probably won't. In those cases it's 
best to pre-slice that decal into separate, convex parts before drawing.

# Future Plans

- Support for text drawing, currently it's not supported, because I'm quite
  lazy.
- Support for concave shapes and holes, these algorithms are complex, so no
  promises

# Building the Example

An example program is included in the repository. This can be built either by
building it manually as a Pixel Game Engine application, or by using the
included meson file.

For the latter, first make sure meson is installed, then run the following
commands:

```sh
meson build
meson compile -C build
./build/olcpgex_viewport
```

In the example, the drawn decal can be moved using the WASD buttons, and rotated
using Q and E
