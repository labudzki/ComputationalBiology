
#ifndef SPATIALSSA_RENDERING_HPP
#define SPATIALSSA_RENDERING_HPP

#include "SpatialSSA_model.hpp"

#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>

void render(Model* model_ptr, bool projection, bool rotate);

#endif