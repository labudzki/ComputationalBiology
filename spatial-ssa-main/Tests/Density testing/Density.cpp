#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>

#include <iostream>
#include <vector>
#include <cmath>

using std::cout;
using std::cin;
using std::endl;

using std::vector;

#include <random>

using std::random_device;
using std::mt19937;
using std::uniform_real_distribution;

// Random number generation
// Generation
random_device rd;  							// Will be used to obtain a seed for the random number engine. As rd() not guaranteed to follow any distribution.
mt19937 gen(rd());							// Standard mersenne_twister_engine seeded with rd(). Is statistically random.
uniform_real_distribution<> dis(0.0, 1.0); 	// Call to "dis(gen)" will now generate a random double in [0,1)


class Circle {
    public:
        double x;
        double y;

        Circle (double x_i, double y_i){
            x = x_i;
            y = y_i;
        }

        double dist(Circle neighbour){
            return sqrt(pow(x-neighbour.x, 2) + pow(y-neighbour.y, 2));
        }
};


int main(){

    int n_particles = 1500;
    int max_attempts = 50;
    bool collision;

    double radius = 10;
    double force_range = 10;

    Circle new_circle(500, 500);
    vector<Circle> circles;

    circles.push_back(new_circle);
    
    for(int i = 0; i<n_particles; i++){

        // cout << "Circle " << i << endl;

        for(int i_attempt = 0; i_attempt < max_attempts; i_attempt++){

            // cout << "Attempt " << i_attempt << endl;

            new_circle.x = 1000*dis(gen);
            new_circle.y = 1000*dis(gen);

            collision = false;

            for(auto circle : circles){
                // cout << circle.x << "   " << circle.y << endl;
                if (circle.dist(new_circle) < 2*radius){
                    collision = true;
                    break;
                }
            }

            if(!collision){
                circles.emplace_back(new_circle.x, new_circle.y);
                break;
            }
        }
    }

    cout << "Number of drawn circles: " << circles.size() << endl;
    cin.get();


    // Rendering

    int width = 1000;
    int height = 1000;

    // double grid_size = 100;

    sf::RenderWindow window(sf::VideoMode(width, height), "Density");

    sf::CircleShape circle(radius);
    sf::RectangleShape grid_line_h(sf::Vector2f(1000, 1));
    sf::RectangleShape grid_line_v(sf::Vector2f(1, 1000));

    while (window.isOpen() ){

        window.clear();

        for(unsigned int i_circle = 0; i_circle < circles.size(); i_circle++){

            circle.setPosition(circles[i_circle].x - radius, circles[i_circle].y - radius);

            if(i_circle == 0){
                circle.setFillColor(sf::Color::Blue);
            }

            else if (circles[i_circle].dist(circles[0]) < 2*radius + force_range){
                circle.setFillColor(sf::Color::Green);
            }

            else {
                circle.setFillColor(sf::Color::White);
            }

            window.draw(circle);

        }

        

        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.display();
    }

    return 0;

}
