#include "SpatialSSA_model.hpp"

#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>


vector<vector<double>> create_colour_scale(int n_colours, int max_n_species){

    vector<vector<double>> colour_scales(n_colours, vector<double>(max_n_species + 1, 0));

    double x_mid, width;
    double x;

    for(int i_colour = 0; i_colour < n_colours; i_colour++){

        x_mid = i_colour * 1.0/(n_colours-1);
        width = 1.0/(n_colours-1);

        for(int i = 0; i < max_n_species+1; i++){

            // x = i/double(max_n_species);
            x = log(i+1)/log(max_n_species+1);
            
            if ((x < x_mid - width) || (x > x_mid + width)){
                colour_scales[i_colour][i] = 0;
            } else {
                colour_scales[i_colour][i] = 0.5 * (cos(M_PI  * (x-x_mid)/(width)) + 1);
            }
        }
    } 

    return colour_scales;
}

vector<vector<double>> create_colour_map(int max_n_species, vector<double> colour_0, vector<double> colour_1, vector<double> colour_2){

    vector<vector<double>> colour_scale = create_colour_scale(3, max_n_species + 1);
    vector<vector<double>> colour(3, vector<double>(max_n_species + 1, 0));

    for(int i_number = 0; i_number < max_n_species+1; i_number++){
        for(int i_channel = 0; i_channel < 3; i_channel++){
            colour[i_channel][i_number] = colour_scale[0][i_number]*colour_0[i_channel]
                                        + colour_scale[1][i_number]*colour_1[i_channel]
                                        + colour_scale[2][i_number]*colour_2[i_channel];
        }
    }

    return colour;
}

vector<vector<double>> create_colour_map(int max_n_species, vector<double> colour_0, vector<double> colour_1, vector<double> colour_2,
                                         vector<double> colour_3, vector<double> colour_4, vector<double> colour_5, vector<double> colour_6,
                                         vector<double> colour_7, vector<double> colour_8, vector<double> colour_9){

    vector<vector<double>> colour_scale = create_colour_scale(10, max_n_species + 1);
    vector<vector<double>> colour(3, vector<double>(max_n_species + 1, 0));

    for(int i_number = 0; i_number < max_n_species+1; i_number++){
        for(int i_channel = 0; i_channel < 3; i_channel++){
            colour[i_channel][i_number] = colour_scale[0][i_number]*colour_0[i_channel]
                                        + colour_scale[1][i_number]*colour_1[i_channel]
                                        + colour_scale[2][i_number]*colour_2[i_channel]
                                        + colour_scale[3][i_number]*colour_3[i_channel]
                                        + colour_scale[4][i_number]*colour_4[i_channel]
                                        + colour_scale[5][i_number]*colour_5[i_channel]
                                        + colour_scale[6][i_number]*colour_6[i_channel]
                                        + colour_scale[7][i_number]*colour_7[i_channel]
                                        + colour_scale[8][i_number]*colour_8[i_channel]
                                        + colour_scale[9][i_number]*colour_9[i_channel];
        }
    }

    return colour;
}


void render(Model* model_ptr, bool projection, bool rotate){

    // In creating grid, store min and max x and y values. Use these to 
    // determine x scale and y scale, use smallest scale to draw.

    // Internal coordinates are in pixels. Starting at top left.
    // Attachment point of shapes (rectangle, circle) is at top left.

    // TODO: in drawing grid, save cell positions, only change colour.

    Reaction* reaction_ptr;
    Cell *cell_ptr, *neighbour_ptr;

    double x0, y0, x_n, y_n, phi, theta, phi_n, theta_n;
    int plot_species = 0;
    plot_species = 4; // Cdc42-GTP MEM
    string LocalComposition;

    // For time printing
    sf::Clock clock;
    sf::Time Current;
    // For cycling through species
    sf::Time LeftClick = clock.getElapsedTime(), RightClick = clock.getElapsedTime();  

    // Determining size of drawn images
    int width = 1600;
    int height = 900;

    sf::RenderWindow window(sf::VideoMode(width, height), "Reaction-Diffusion");

    sf::RectangleShape bg(sf::Vector2f(width, height));
    //bg.setSize(sf::Vector2f(100, 50));
    bg.setPosition(0, 0);
    bg.setFillColor(sf::Color(255, 255, 255));

    // Area of drawable canvas does not include the boder
    double border = 50.f, s_border = 2.5;
    
    // For printing text;
    int text_size = 18;

    // Idea: the cells in the model have certain coordinates. We use these coordinates to draw a map of the network. 
    // We determine the size scaling by comparing the maximum distance in internal coordinates in the x and y direction,
    // with the max size of the window (minus a border size).
    
    double v_scale, h_scale, scale, s_size, l_size;
    width -= 2*border;
    height -= 2*border;

    // Total distance in internal coordinates without shape size incorporated.
    double tot_x, tot_y;
    double min_x = model_ptr->cells[0].pos.x, max_x = model_ptr->cells[0].pos.x;
    double min_y = model_ptr->cells[0].pos.y, max_y = model_ptr->cells[0].pos.y;
    double dist, min_dist;

    int max_n_species = 0;

    // Minimal distance per dimension
    int draw_dim = 2;

    for(int i_cell = 0; i_cell < model_ptr->n_cells; i_cell++) {

        cell_ptr = &(model_ptr->cells[i_cell]);
        if(cell_ptr->dim != draw_dim){
            continue;
        }

        // Finding capacity to determine colour scale
        if (int(cell_ptr->capacity) + 1 > max_n_species){
            max_n_species = int(cell_ptr->capacity) + 1;
        }
        
        // Finding minimal and maximal real or projection coordinates (so not in pixels yet)
        if (projection) {

            x0 = atan2(cell_ptr->pos.y, cell_ptr->pos.x);
            y0 = acos(cell_ptr->pos.z/sqrt(pow(cell_ptr->pos.x, 2) + pow(cell_ptr->pos.y, 2) + pow(cell_ptr->pos.z, 2)));
        } else {
            x0 = cell_ptr->pos.x;
            y0 = cell_ptr->pos.y;
        }

        min_x = min(min_x, x0);
        max_x = max(max_x, x0);

        min_y = min(min_y, y0);
        max_y = max(max_y, y0);


        // Finding minimal distance in real coordinates
        for(unsigned int i_neighbour = 0; i_neighbour < cell_ptr->n_neighbours; i_neighbour++){

            // cout << "Index neighbour: " << i_cell << endl;
            neighbour_ptr = &model_ptr->cells[cell_ptr->neighbours[i_neighbour]];

            if(neighbour_ptr->dim != draw_dim){
                continue;
            }

            if (projection) {
                x_n = atan2(neighbour_ptr->pos.y, neighbour_ptr->pos.x);
                y_n = acos(neighbour_ptr->pos.z/sqrt(pow(neighbour_ptr->pos.x, 2) + pow(neighbour_ptr->pos.y, 2) + pow(neighbour_ptr->pos.z, 2)));
            } else{
                x_n = neighbour_ptr->pos.x;
                y_n = neighbour_ptr->pos.y;
            }

            dist = sqrt(pow(x0 - x_n, 2) + pow(y0 - y_n, 2));

            if (i_cell == 0 && i_neighbour == 0){
                min_dist = dist;
            }
            min_dist = min(min_dist, dist);
        }
    }

    max_n_species = 3000;
    // max_n_species = 200;
    // max_n_species = 50;
    // max_n_species = 20;

    // cout << "Distance x: " << min_x << " to " << max_x << endl;
    // cout << "Distance y: " << min_y << " to " << max_y << endl;
    // cout << "Shortest distance: " << min_dist << endl;

    tot_x = (max_x - min_x) + min_dist;
    tot_y = (max_y - min_y) + min_dist;

    // cout << tot_x << "   " << tot_y << endl;
    // 2pi, pi

    // Use height, not width, because I want part of the screen for other stuff
    v_scale = height/tot_x;
    h_scale = height/tot_y;
    scale = min(v_scale, h_scale);

    // The maximum size individual cells can be drawn at, without them touching.
    s_size = scale*min_dist - 2*s_border;
    s_size = max(scale*min_dist - 2*s_border, 8.);

    double s_edge = 1;
    // cout << "Min dist: " << min_dist << endl;  
    // cout << "Total x: " << tot_x << endl;
    // cout << "Total y: " << tot_y << endl;
    // cout << "Scale: " << scale << endl;
    // cout << "Shape size: " << s_size << endl;
    // cin.get();
 
    // Legenda variables.
    double count_species, factor;

    int N_legend = 15;

    // The size of the legenda boxes.
    l_size = double(height)/(N_legend) - 2*s_border;
    sf::RectangleShape legend_box(sf::Vector2f(l_size, l_size));

    // Color scales
    // Viridis colour map
    vector<double> colour_0{ 68, 13, 86 };
    vector<double> colour_1{ 32, 143, 141 };
    vector<double> colour_2{ 247, 229, 31 };
    
    // Plasma colour map
    // vector<double> colour_0{ 16, 7, 135 };
    // vector<double> colour_1{ 195, 62, 127 };
    // vector<double> colour_2{ 239, 247, 34 };

    // Inferno colour map
    // vector<double> colour_0{ 0, 0, 4 };
    // vector<double> colour_1{ 178, 50, 89 };
    // vector<double> colour_2{ 250, 250, 160 };

    // Magma colour map
    // vector<double> colour_0{ 0, 0, 6 };
    // vector<double> colour_1{ 172, 51, 123 };
    // vector<double> colour_2{ 250, 250, 190 };

    vector<vector<double>> colour = create_colour_map(max_n_species, colour_0, colour_1, colour_2);

    // custom colour map
    // vector<double> colour_9{ 195, 206, 118 };
    // vector<double> colour_8{ 143, 183, 132 };
    // vector<double> colour_7{ 90, 158, 148 };
    // vector<double> colour_6{ 36, 132, 154 };
    // vector<double> colour_5{ 45, 116, 174 };
    // vector<double> colour_4{ 53, 97, 184 };
    // vector<double> colour_3{ 63, 79, 195 };
    // vector<double> colour_2{ 96, 71, 173 };
    // vector<double> colour_1{ 129, 64, 152 };
    // vector<double> colour_0{ 164, 58, 130 };

    // vector<vector<double>> colour = create_colour_map(max_n_species, colour_0, colour_1, colour_2, colour_3, colour_4,
    //                                                   colour_5, colour_6, colour_7, colour_8, colour_9);


    string state_value, messages;
    
    // Setting known positions of the cells
    vector<sf::RectangleShape> cell_shapes(model_ptr->n_cells, sf::RectangleShape(sf::Vector2f(s_size, s_size)));
    vector<sf::RectangleShape> edges(model_ptr->n_cells * MAX_NEIGHBOURS, sf::RectangleShape(sf::Vector2f(1, 1)));
    sf::Vector2f rectanglePosition; 
    sf::Vector2i localPosition;
    int n_edges = 0, neighbour_idx;


    // Determining the pixel coordinates of all the shapes and storing them.
    for(int i_cell = 0; i_cell < model_ptr->n_cells; i_cell++) {

        cell_ptr = &(model_ptr->cells[i_cell]);
        
        if (projection) {
            phi = atan2(cell_ptr->pos.y, cell_ptr->pos.x);
            theta = acos(cell_ptr->pos.z/sqrt(pow(cell_ptr->pos.x, 2) + pow(cell_ptr->pos.y, 2) + pow(cell_ptr->pos.z, 2)));

            x0 = border + s_border + scale*(phi - min_x);
            y0 = border + s_border + 2*scale*(theta - min_y);
        } else {
            // Top right coordinate of the current cell.
            x0 = border + s_border + scale*(cell_ptr->pos.x - min_x);
            y0 = border + s_border + scale*(cell_ptr->pos.y - min_y);
        }
        
        cell_shapes[i_cell].setPosition(x0,y0);
    }

    // The actual shapes of the cells
    vector<sf::ConvexShape> hexagons(model_ptr->n_cells, sf::ConvexShape());
    if(projection){
        for (int i_cell = 0; i_cell < model_ptr->n_cells; i_cell++) {
            cell_ptr = &model_ptr->cells[i_cell];

            hexagons[i_cell].setPointCount(cell_ptr->n_neighbours);

            phi_n = atan2(cell_ptr->pos.y, cell_ptr->pos.x);
            theta_n = acos(cell_ptr->pos.z/sqrt(pow(cell_ptr->pos.x, 2) + pow(cell_ptr->pos.y, 2) + pow(cell_ptr->pos.z, 2)));

            for (unsigned int i_point = 0; i_point < cell_ptr->n_neighbours; i_point++){

                phi = atan2(cell_ptr->vertices[i_point].y, cell_ptr->vertices[i_point].x);
                theta = acos(cell_ptr->vertices[i_point].z/sqrt(pow(cell_ptr->vertices[i_point].x, 2) + pow(cell_ptr->vertices[i_point].y, 2) + pow(cell_ptr->vertices[i_point].z, 2)));

                if(phi - phi_n > M_PI){phi -= 2*M_PI;}
                if(phi - phi_n < -M_PI){phi += 2*M_PI;}
                x0 = scale*(phi - phi_n);
                // cout << x0 << " ";
                y0 = 2*scale*(theta - theta_n);

                hexagons[i_cell].setPoint(i_point, sf::Vector2f(x0, y0));
            }

            // cin.get();
            
            // Scaling hexagons, so that you have a bit of a border
            // hexagons[i_cell].setScale(0.8, 0.8);
            hexagons[i_cell].move(cell_shapes[i_cell].getPosition() + sf::Vector2f(s_size/2, s_size/2));
            hexagons[i_cell].scale(.8, .8);
        }
    }

    // Determining where the edges go for the spherical grid layout.
    if(projection){
        for (int i_cell = 0; i_cell < model_ptr->n_cells; i_cell++) {

            cell_ptr = &(model_ptr->cells[i_cell]);
            if(cell_ptr->dim != draw_dim){
                continue;
            }

            for (unsigned int i_neighbour = 0; i_neighbour < cell_ptr->n_neighbours; i_neighbour++) {
                neighbour_idx = cell_ptr->neighbours[i_neighbour];

                rectanglePosition = cell_shapes[i_cell].getPosition();
                phi = rectanglePosition.x, theta = rectanglePosition.y;

                rectanglePosition = cell_shapes[neighbour_idx].getPosition();
                phi_n = rectanglePosition.x, theta_n = rectanglePosition.y;

                if(sqrt(pow(phi_n - phi, 2) + pow(theta_n - theta, 2)) > 2./3*height) {continue;}

                edges[n_edges].setPosition(phi + s_size/2, theta + s_size/2);
                edges[n_edges].setSize(sf::Vector2f(sqrt(pow(phi_n - phi, 2) + pow(theta_n - theta, 2)), s_edge));
                edges[n_edges].setRotation(180/M_PI * atan2(theta_n - theta, phi_n - phi));

                edges[n_edges].setFillColor(sf::Color(100, 100, 100));

                n_edges++;
            }   
        }
    }

    
    sf::Text text;
    text.setCharacterSize(text_size);
    text.setFillColor(sf::Color::White);

    sf::Font font;
    if (!font.loadFromFile("Fonts/arial.ttf")) {
        cout << "Could not load font..." << endl;
    }

    text.setFont(font);
	text.setOutlineColor(sf::Color::Black);	        
	text.setOutlineThickness(2);

    // Variables for turning
    double pos_x, pos_z;
    bool turn = false;
    int d_tau = 0, tau = 20;
    double omega = 0.05*2*M_PI*tau/1000;

	
    while (window.isOpen() ){

        window.clear();

        // draw background
        //window.draw(bg);

        // Drawing text

        messages = "Elapsed clock time: " + to_string(clock.getElapsedTime().asSeconds()) + "\n";   // Elapsed time
        messages += "Elapsed simulation time: " + to_string(model_ptr->t) + "\n";                    // Elapsed simulation time
        messages += "Iteration number: " + to_string(model_ptr->i_iter) + "\n\n";                      // Elapsed iterations

        // Counts of all types of reactions
        for(int i_reaction = 0; i_reaction < model_ptr->n_reactions; i_reaction++){
            reaction_ptr = &(model_ptr->reactions[i_reaction]);
            messages += reaction_ptr->reaction_name + " count: " + to_string(reaction_ptr->count) + "\n";
        }
        messages += "\n";

        // Counts of all types of reactions
        for(int i_reaction = 0; i_reaction < model_ptr->n_boundary_reactions; i_reaction++){
            reaction_ptr = &(model_ptr->boundary_reactions[i_reaction]);
            messages += reaction_ptr->reaction_name + " count: " + to_string(reaction_ptr->count) + "\n";
        }
        messages += "\n";

        // Counts of all types of diffusion
        for(int i_species = 0; i_species < model_ptr->n_species; i_species++){
            reaction_ptr = &(model_ptr->diffusion_reactions[i_species]);
            messages += reaction_ptr->reaction_name + " count: " + to_string(reaction_ptr->count) + "\n";
        }

        messages += "\n";
        messages += "Showing species " + model_ptr->species[plot_species%model_ptr->n_species] + "\n";

        // Drawing messages
        text.setString(messages);
        text.setPosition(border + height + 0.5*border + (l_size + 2*s_border) + 0.5*border, border);
        window.draw(text);

        x0 = border + height + 0.5*border; 

		for (int i_legend = 0; i_legend < N_legend; i_legend++) {
		
			y0 = border + s_border + i_legend * (l_size + 2*s_border); 
		
	    	legend_box.setPosition(x0,y0);
	    	text.setPosition(x0,y0);
	    	
	    	factor = pow(max_n_species,i_legend / double(N_legend-1));
	    	text.setString(to_string(int(round(factor))));

            legend_box.setFillColor(sf::Color(colour[0][int(factor)], colour[1][int(factor)], colour[2][int(factor)]));
			window.draw(legend_box);
			window.draw(text);
		}

        // Getting position of mouse cursor.
		localPosition = sf::Mouse::getPosition(window);

        // Drawing edges
        for(int i_edge = 0; i_edge < n_edges; i_edge++) {
            window.draw(edges[i_edge]);
        }

        
        // Is it time to turn?
        if(rotate){
            if(turn){turn = false;}
            if(d_tau < clock.getElapsedTime().asMilliseconds()/tau){   
                d_tau++;
                turn = true;
            }
        }

        // Drawing cells
        for(int i_cell = 0; i_cell < model_ptr->n_cells; i_cell++) {

            cell_ptr = &(model_ptr->cells[i_cell]);
            if(cell_ptr->dim != draw_dim){
                continue;
            }

            // Rotating cells
            pos_x = cell_ptr->pos.x; 
            pos_z = cell_ptr->pos.z;

            if(turn && rotate){
                x0 = border + s_border + scale*(cos(omega*d_tau)*pos_x + sin(omega*d_tau) * pos_z - min_x);
                y0 = border + s_border + scale*(cell_ptr->pos.y - min_y);
        
                cell_shapes[i_cell].setPosition(x0,y0);
            }

            rectanglePosition = cell_shapes[i_cell].getPosition();

            count_species = cell_ptr->state[plot_species%model_ptr->n_species];
            // count_species = cell_ptr->idx;
            cell_shapes[i_cell].setFillColor(sf::Color(colour[0][int(count_species)], colour[1][int(count_species)], colour[2][int(count_species)]));
            
            if(!rotate || -sin(omega*d_tau)*pos_x + cos(omega*d_tau) * pos_z > 0){
                window.draw(cell_shapes[i_cell]);
            }

            // hexagons[i_cell].setFillColor(sf::Color(colour[0][int(count_species)], colour[1][int(count_species)], colour[2][int(count_species)]));
            // window.draw(hexagons[i_cell]);


            // If mouse in the current cell, display current contents of that cell.
            if (localPosition.x > rectanglePosition.x and localPosition.x < (rectanglePosition.x + s_size)
                and localPosition.y > rectanglePosition.y and localPosition.y < (rectanglePosition.y + s_size)) {
                LocalComposition = "Cell " + to_string(cell_ptr->idx) + "\n";
                if(model_ptr->verbose){
                    LocalComposition += "Neighbours " + to_string(cell_ptr->n_neighbours) + "\n";
                    LocalComposition += "X " + to_string(cell_ptr->pos.x) + "\n";
                    LocalComposition += "Y " + to_string(cell_ptr->pos.y) + "\n";
                    LocalComposition += "Z " + to_string(cell_ptr->pos.z) + "\n";
                }
                for(int i_species = 0; i_species < model_ptr->n_species; i_species++){
	                LocalComposition += model_ptr->species[i_species] + " " + to_string(cell_ptr->state[i_species]) + "\n";
	            }

				text.setString(LocalComposition);
				text.setPosition(localPosition.x+30,localPosition.y-30);        
	        }

        }

        window.draw(text);

        // For switching which species is currently rendered.
        // It checks whether the left mouse button has been clicked. If it is clicked, and time has advanced more than 200 ms,
        // then make change to the plotted species and reset click time. Ensures held down click does not advance species.
		if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
			Current = clock.getElapsedTime(); 
			if (Current.asMilliseconds() - LeftClick.asMilliseconds() > 200) {
				plot_species++;
				LeftClick = clock.getElapsedTime(); 
			}
		}
		if (sf::Mouse::isButtonPressed(sf::Mouse::Right)) {
			Current = clock.getElapsedTime(); 
			if (Current.asMilliseconds() - RightClick.asMilliseconds() > 200) {
                if(plot_species == 0) {plot_species = model_ptr->n_species;}
				plot_species--;
				RightClick = clock.getElapsedTime(); 
			}
		}
    	
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.display();
    }
}