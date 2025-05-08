#include "SpatialSSA_model.hpp"

// // For (trying to set) a random seed
// #include <cstdlib>
// #include <ctime>
// #include <cstdio>

#include <random>

using std::random_device;
using std::mt19937;
using std::uniform_real_distribution;

// Timing of code

#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::duration_cast;

high_resolution_clock timing_clock;
auto t1 = high_resolution_clock::now();
auto t2 = high_resolution_clock::now();
auto t_start = high_resolution_clock::now();
auto t_now = high_resolution_clock::now();
// auto ta = high_resolution_clock::now();
// auto tb = high_resolution_clock::now();

// Getting number of milliseconds as an integer.
auto total_time = seconds(0);
auto affected_cell_time = milliseconds(0);
auto neighbour_cell_time = milliseconds(0);
auto find_mu_time = milliseconds(0);
auto update_time = milliseconds(0);
auto record_time = milliseconds(0);

// Random number generation

// For compilers which have a random device, we can simply use that to set the seed
// Will be used to obtain a seed for the random number engine. As rd() not guaranteed to follow any distribution.
random_device rd;  							
// Standard mersenne_twister_engine seeded with rd(). Is statistically random.
mt19937 gen(rd());

// For compilers which do not have an intrinsic source of randomness, set seed as time
// mt19937 gen(static_cast<long unsigned int>(time(0)));							
uniform_real_distribution<> dis(0.0, 1.0); 	// Call to "dis(gen)" will now generate a random double in [0,1)

void tell(string message, bool fail){
    cout << message << endl;
    cin.get();
    if(fail){
        exit(1);
    }
}

Reaction::Reaction() {

    n_reaction_species = 0;
    n_reactants = 0;

    count = 0;
}

Cell::Cell(){

    vertices.reserve(MAX_NEIGHBOURS);
    neighbours.resize(MAX_NEIGHBOURS);
    affected_events.resize(MAX_REACTIONS + 2*MAX_NEIGHBOURS*MAX_NEIGHBOURS*MAX_SPECIES);
    // Reactions in cell, diffusions with neighbour cells, and diffusions of neighbours

    size = 0;

    n_neighbours = 0;
    n_affected_events = 0;

    for(int i = 0; i < MAX_SPECIES; i++){
        state[i] = 0;
    }
}

double Cell::distance(const Cell &neighbour){

    return pos.distance(neighbour.pos);
}


void Cell::set_species_number(int species_idx, int value){

    // Add check that species exists

    state[species_idx] = value;
}


void Cell::add_neighbour(int neighbour_idx, double distance, double interface){
    // Adding neighbour of the same dimensionality type

    if(n_neighbours == MAX_NEIGHBOURS){
        cout << "Number of neighbours of cell " << idx << " (" << n_neighbours + 1 << ") exceeds maximum number of neighbours" << endl;
        cin.get();
        exit(1);
    }
    neighbours[n_neighbours] = neighbour_idx;
    neighbour_dist[neighbour_idx] = distance;
    neighbour_interface[neighbour_idx] = interface;
    n_neighbours++;

}

void Cell::add_neighbour_alt_dim(int neighbour_idx, double interface){
    // Adding neighbour of a different dimensionality type
    
    if(n_neighbours == neighbours.size()){
        cout << "Number of neighbours of cell " << idx << " (" << n_neighbours + 1 << ") exceeds maximum number of neighbours" << endl;
        cin.get();
        exit(1);
    }

    neighbours[n_neighbours] = neighbour_idx;
    neighbour_dist[neighbour_idx] = 0;
    neighbour_interface[neighbour_idx] = interface;
    n_neighbours++;

}

bool Cell::has_neighbour(int neighbour_idx){

    for(unsigned int i_neighbour = 0; i_neighbour < n_neighbours; i_neighbour++){
        if(neighbours[i_neighbour] == neighbour_idx){
            return true;
        }
    }
    return false;

}

void Cell::add_affected_event(Event* event_ptr){

    if(n_affected_events >= affected_events.size()){
        cout << "Number of events exceeds maximum number of affected events (" << affected_events.size() << ")" << endl;
        cin.get();
        exit(1);
    }
    affected_events[n_affected_events] = event_ptr;
    n_affected_events++;
}

Model::Model(){

    phase_separate = false;
    crowding       = false;


    // Phase seperation parameters
    particles_per_area = 40000.;
    particles_per_volume = 8000000.;
    chi = 0; kappa = 0;

    // Crowding parameters
    free_area_fraction = 1;
    free_volume_fraction = 1;

    verbose = false;
    timing = false;
    initial_recording = true;

    record_indv_cells = record_indv_species = record_reactions = true;
    record_propensities = true;

    n_repeats = 1;

    t_max = 0;
    n_iterations = 0;

    i_run = 0;

    n_cells = 0;
    n_reactions = 0;
    n_species = 0;
    n_events = 0;

    t = 0;

    for(int i_event = 0; i_event < MAX_EVENTS; i_event++){
        propensities[i_event] = 0;
    }

    for(int i_species = 0; i_species < MAX_SPECIES; i_species++){
        initial_species_number[i_species] = 0;
    }

    // Setting project date. Time the model is run
    time_t start_time = time(0);
    tm *ltm = localtime(&start_time);

    project_date = to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday) + "_"
                 + to_string(ltm->tm_hour) + "-" + to_string(ltm->tm_min);
}

void Model::reset(){

    project_name = "./" + dir_name + "/" + project_date + "_run_" + to_string(i_run) + "_" + run_name;

    initial_recording = true;

    n_events = 0;
    t = 0;
    i_iter = 0;

    // Set the free_size fraction per dimension.
    free_size_fraction[2] = free_area_fraction;
    free_size_fraction[3] = free_volume_fraction;

    for(int i_reaction = 0; i_reaction < n_reactions; i_reaction++){
        reactions[i_reaction].count = 0;
    }

    for(int i_dif_reaction = 0; i_dif_reaction < n_species; i_dif_reaction++){
        diffusion_reactions[i_dif_reaction].count = 0;
    }

    for(int i_bnd_reaction = 0; i_bnd_reaction < n_boundary_reactions; i_bnd_reaction++){
        boundary_reactions[i_bnd_reaction].count = 0;
    }

    // Reset affected events
    for(int i_cell = 0; i_cell < n_cells; i_cell++){
        cells[i_cell].n_affected_events = 0;
    }

    // Redistribute species
    for(int i_species = 0; i_species < n_species; i_species++){
        distribute_particles(initial_species_number[i_species], i_species);
    }


    F_max = 0;
    n_reject = 0;
    n_recalculate = 0;

    mu_history.clear();
    mu_history.reserve(break_iter);

}


void Model::make_particle_spot(const int cell_idx, const double spot_radius, const int count, const int i_species){
    
    int i_cell;
    vector<int> cells_to_visit;
    set<int> visited_cells;
    set<int>::iterator itr;

    Cell *center_cell_ptr = &cells[cell_idx];
    Cell *cur_cell_ptr, *neighbour_ptr;

    cells_to_visit.push_back(cell_idx);
    visited_cells.insert(cell_idx);
    
    while(cells_to_visit.size() != 0){
        
        // Get last element and remove it from list
        i_cell = cells_to_visit.back();
        cells_to_visit.pop_back();
            
        // cout << "Looking at cell " << i_cell << endl;
        // cout << "Number of cells to visit: " << cells_to_visit.size() << endl;

        cur_cell_ptr = &cells[i_cell];

        // cin.get();

        for (unsigned int i_neighbour = 0; i_neighbour < cur_cell_ptr->n_neighbours; i_neighbour++){
            neighbour_ptr = &cells[cur_cell_ptr->neighbours[i_neighbour]];

            // cout << "Neighbour cell " << neighbour_ptr->idx << endl;
            // cout << "Distance from center: " << center_cell_ptr->distance(*neighbour_ptr) << endl;

            // If cell is close enough and not yet visited
            if ((neighbour_ptr->dim == center_cell_ptr->dim) && (center_cell_ptr->distance(*neighbour_ptr) < spot_radius)
                 && (visited_cells.find(neighbour_ptr->idx) == visited_cells.end())){
                    
                    // cout << "Added cell " << neighbour_ptr->idx << endl;
                    cells_to_visit.push_back(neighbour_ptr->idx);
                    visited_cells.insert(neighbour_ptr->idx);
            }
        }

        // cin.get();
    }


    itr = visited_cells.begin();

    for(int i_spot_particle = 0; i_spot_particle < count; i_spot_particle++){

        i_cell = *itr;

        // cout << "Adding particle to cell " << i_cell << endl;
        // cin.get();

        cells[i_cell].state[i_species] += 1;
        cells[i_cell].content += 1;
        cells[i_cell].phi = (cells[i_cell].content)/cells[54].capacity;

        itr++;
        if (itr == visited_cells.end()){itr = visited_cells.begin();}
    }
}

void Model::distribute_particles(const int count, const int i_species){

    double p_cell, r;
    double tot_capacity = 0;
    int i_cell, species_dim = species_dimensions[i_species];

    // cout << "Species name: " << species[i_species] << endl;
    // cout << "Species dim: " <<species_dim << endl;
    // cin.get();

    if(n_cells_per_dim[species_dim] == 0){
        tell("Number of cells in model must be greater than 0", 1);
    }

    // Clear cells of said species
    for(int j_cell = 0; j_cell < n_cells; j_cell++){
        cells[j_cell].state[i_species] = 0;
    }
    
    // Test to see if there is enough space
    // Todo: change it so that it sums over all cells and sums up space.
    // It must not sum up space, but still available space.
    for(int j_cell = 0; j_cell < n_cells; j_cell++){
        if (cells[j_cell].dim == species_dim){
            // cout << "Capacity per cell: " << cells[j_cell].capacity << endl;
            tot_capacity += cells[j_cell].capacity;
            
        }
    }

    if(count >= tot_capacity){
        cout << "Number of particles to be added exceeds number of available spaces." << endl;
        cout << "Total capacity: " << tot_capacity << endl;
        cout << "Particles to be added: " << count << endl;
        cin.get();
        exit(1);
    }


    // Add particles randomly
    p_cell = 1/double(n_cells_per_dim[species_dim]);
    double p_cum;
    bool added_particle;

    for(int i_particle = 0; i_particle < count; i_particle++){
        // cout << "Adding particle " << i_particle << endl;

        added_particle = false;
        
        
        // Loop so that particles only get added if there is space left
        while(!added_particle){

            r = dis(gen);
            p_cum = 0;

            // Loop over cells, and return cell index with correct species based on random number
            for (i_cell = 0; i_cell < n_cells; i_cell++){
                if (cells[i_cell].dim == species_dim){
                    p_cum += p_cell;
                    if (p_cum > r){break;}
                }
            }

            // Check if there is space. If not, try again
            if(cells[i_cell].content <= 0.9 * cells[i_cell].capacity){
                cells[i_cell].state[i_species] += 1;
                cells[i_cell].content += 1;
                cells[i_cell].phi = (cells[i_cell].content)/cells[i_cell].capacity;
                added_particle = true;
            }
        }
    }
}


void Model::add_rate_constant(const string rate_constant_name){

    if(n_rate_constants == MAX_REACTIONS+MAX_SPECIES){
        tell("Number of rate constants exceeds maximum number of reactions (and hence constants)", 1);
    }

    if(rate_constant_map.count(rate_constant_name)){
        cout << "Rate constant already exists";
        cin.get();
        exit(1);
    }

    rate_constant_map[rate_constant_name] = n_rate_constants;
    rate_constant_names[n_rate_constants] = rate_constant_name;
    n_rate_constants++;

    // Logging
    cout << "Added rate constant " << rate_constant_names[n_rate_constants - 1] << " with index " << rate_constant_map[rate_constant_names[n_rate_constants - 1]] << endl;;

    // // Check if rate constant exists already
    // for(int i = 0; i < n_rate_constants; i++){
    //     if(rate_constant_name == rate_constant_names[i]){
    //         cout << "Rate constant already exists";
    //         cin.get();
    //         exit(1);
    //     }
    // }

    // rate_constant_names[n_rate_constants] = rate_constant_name;
    // n_rate_constants++;
}


void Model::add_species(const string species_name, const int dimension, const string dif_coefficient_name){

    Reaction* reaction_ptr;

    int dif_coefficient_idx = get_rate_idx(dif_coefficient_name);

    if(n_species == MAX_SPECIES){
        tell("Number of species exceeds maximum number of species", 1);
    }

    // Check if species exists already
    for(int i = 0; i < n_species; i++){
        if(species_name == species[i]){
            cout << "Species already exists";
            return;
        }
    }

    species[n_species] = species_name;
    species_dimensions[n_species] = dimension;

    reaction_ptr = &diffusion_reactions[n_species];
    reaction_ptr->k_idx = dif_coefficient_idx;

    // In both cells, the same species is involved
    reaction_ptr->reaction_species[0] = n_species;
    reaction_ptr->reaction_species[1] = n_species;
    reaction_ptr->n_reaction_species = 2;

    // In one cell, a species disappears, which appears in the other cell
    reaction_ptr->state_change[0] = -1;
    reaction_ptr->state_change[1] = 1;

    // Add reactants
    reaction_ptr->n_reactants = 1;
    reaction_ptr->reactants[0] = 0;

    reaction_ptr->reaction_name = species_name + " diffusion";
    reaction_ptr->is_diffusion_like = true;

    n_species++;

    // Logging
    cout << "Added species " << species[n_species - 1] << " with dimension " << species_dimensions[n_species - 1] << ", diffusion constant index " << reaction_ptr->k_idx << endl;;
    cout << endl;
}

void Model::add_reaction(const string reaction_species[], const int state_changes[], const string rate_constant_name, const int n_reaction_species, const string reaction_name){

    if(n_reactions == MAX_REACTIONS){
        tell("Number of reactions exceeds maximum number of reactions", 1);
    }

    if(n_reaction_species > MAX_REACTION_SIZE){
        cout << "Number of reactants in reaction " << reaction_name << " exceeds maximum number of reactants" << endl;
        cin.get();
        exit(1);
    }


    // Temporary variables
    Reaction* reaction_ptr;
    bool is_boundary_reaction = false;
    int species_idx, species_dim, rate_constant_idx;

    // Determine whether the reaction is normal, or that it crosses dimensionality boundaries
    species_dim = species_dimensions[get_species_idx(reaction_species[0])];
    for(int i = 0; i<n_reaction_species; i++){
        // If the previous species dimension is not equal
        if (species_dim != species_dimensions[get_species_idx(reaction_species[i])]){
            is_boundary_reaction = true;
            break;
        }
        species_dim = species_dimensions[get_species_idx(reaction_species[i])];
    }

    // Get pointer to current reaction
    if(is_boundary_reaction){
        reaction_ptr = &boundary_reactions[n_boundary_reactions];
        n_boundary_reactions++;
    } else{
        reaction_ptr = &reactions[n_reactions];
        n_reactions++;
    }
    
    
    // Loop over all species in the reaction to see if they exist, and if so, add their index to the reaction class
    for(int i = 0; i<n_reaction_species; i++){

        if(state_changes[i] < -1){
            cout << "The state change for species " << reaction_species[i] << " in reaction " << reaction_name << " is less than -1." << endl;
            cout << "This model supports linear terms as is, and quadratic terms (in the same reaction species) by commenting in an additional part in the get_event_propensity function" << endl;
            cout << "If the term is more than quadratic, you have to add in support for this manually, by adding additional cases in calculate_propensities" << endl;
            cin.get();
            exit(1);
        }

        species_idx = get_species_idx(reaction_species[i]);

        reaction_ptr->reaction_species[i] = species_idx;
        reaction_ptr->state_change[i] = state_changes[i];

        // Determining which species are reactants
        // Todo: add better support for catalysis, by allowing species to have a state change of 0,
        // and still contribute to the propensity.
        if(state_changes[i] < 0){
            reaction_ptr->reactants[reaction_ptr->n_reactants] = i;
            reaction_ptr->n_reactants++;
        } 
    }

    rate_constant_idx = get_rate_idx(rate_constant_name);

    reaction_ptr->k_idx = rate_constant_idx;
    reaction_ptr->n_reaction_species = n_reaction_species;

    reaction_ptr->reaction_name = reaction_name;
    reaction_ptr->is_diffusion_like = false;

    // Logging for safety
    cout << "Added reaction " << reaction_ptr->reaction_name << endl;
    for (int i_reaction_species = 0; i_reaction_species < reaction_ptr->n_reaction_species; i_reaction_species++){
        cout << "Species " << species[reaction_ptr->reaction_species[i_reaction_species]] << " , change: " << reaction_ptr->state_change[i_reaction_species] << endl;
    }

    cout << "Reactants are:" << endl;
    for (int i_reactant = 0; i_reactant < reaction_ptr->n_reactants; i_reactant++){
        cout << "Species " << species[reaction_ptr->reaction_species[reaction_ptr->reactants[i_reactant]]] << endl;
    }

    // Some reactions also must be treated as diffusion in phase seperation, such as membrane detachment
    if(reaction_ptr->n_reaction_species == 2 && reaction_ptr->n_reactants == 1 &&
        species_dimensions[reaction_ptr->reaction_species[reaction_ptr->reactants[0]]] == 2){

        reaction_ptr->is_diffusion_like = true;
        cout << "Added reaction " << reaction_name << " to diffusion-like reactions" << endl;
    }

    cout << endl;

}

int Model::get_species_idx(const string species_name){

    int species_idx = -1;

    for(int i = 0; i < n_species; i++){
            if(species[i] == species_name){
                species_idx = i;
            }
    }

    if(species_idx == -1){
        cout << "Reaction species " << species_name << " was not yet defined" << endl;
        cin.get();
        exit(1);
    }

    return species_idx;
}

int Model::get_rate_idx(const string rate_constant_name){

    if(!rate_constant_map.count(rate_constant_name)){
        cout << "Rate constant " << rate_constant_name << " was not yet defined" << endl;
        cin.get();
        exit(1);
    }

    return rate_constant_map[rate_constant_name];

    // int rate_idx = -1;

    // for(int i = 0; i < n_rate_constants; i++){
    //         if(rate_constant_names[i] == rate_constant_name){
    //             rate_idx = i;
    //         }
    // }

    // if(rate_idx == -1){
    //     cout << "Rate constant " << rate_constant_name << " was not yet defined" << endl;
    //     cin.get();
    //     exit(1);
    // }

    // return rate_idx;
}

void Model::add_species_sweep(const int i_species, const int n_start, const int n_end, const int n_points, const string sweep_type){

    if(n_points <= 1){
        tell("Number of points in the sweep must be bigger than 1", 1);
    }

    vector<int> sweep;

    if(sweep_type == "lin"){
        for(int i_point = 0; i_point < n_points; i_point++){
            sweep.push_back(int(i_point*double(n_end-n_start)/(n_points-1)) + n_start);
        }
    }
    else if(sweep_type == "log"){

        if(n_start == 0){
            tell("In log sweep, start number cannot be 0", 1);
        }

        for(int i_point = 0; i_point < n_points; i_point++){
            sweep.push_back(int(n_start * pow(double(n_end)/n_start, double(i_point)/(n_points-1))));
        }
    }

    // Add to model member
    sweeps_species_number.push_back(sweep);
    sweeps_species_idx.push_back(i_species);

}

// inline double Model::calculate_gradient_energy(){

//     double F = 0;

//     // cout << "Number of faces: " << faces.size() << endl;

//     for (Face face : faces){
//         // Loop over vertices
//         for(int i = 0; i<3; i++){

//             // cout << "Density: " << face.cells[i]->phi << endl;
//             // cout << "This term: " << pow(face.lengths[i]*face.cells[i]->phi, 2) - 2*face.theta_terms[i] * (face.cells[(i+1)%3]->phi) * (face.cells[(i+2)%3]->phi) << endl;
            
//             F += pow(face.lengths[i]*face.cells[i]->phi, 2) - 2*face.theta_terms[i] * (face.cells[(i+1)%3]->phi) * (face.cells[(i+2)%3]->phi);
            
//         }

//         // cout << "Face " << i_face << endl;
//         // cout << "This face: " << dF << endl;
//         // cout << "Total F: " << F << endl;
//         // cin.get();
//     }
//     // Prefactor from averaging over two phi's, squared.
//     // Here, F is only the gradient squared.
//     F /= 4;

//     // Energy scaling factor applied to gradient squared.
//     return kappa/2 * F;
// }

// inline double Model::calculate_gradient_energy_dif(const Event* event_ptr){

//     Cell *cell_ptr_1, *cell_ptr_2;
//     double del_phi_1, del_phi_2;

//     double del_f = 0;   // Total energy density
//     double del_F = 0;   // Total energy

//     cell_ptr_1 = &cells[event_ptr->target_cells[0]];
//     cell_ptr_2 = &cells[event_ptr->target_cells[1]];

//     del_phi_1 = -1/cell_ptr_1->capacity;   // Todo: change in next iteration to be particle dependent and cell size dependent
//     del_phi_2 = 1/cell_ptr_2->capacity;

//     del_F += calculate_gradient_energy_dif_mesh(cell_ptr_1->phi, del_phi_1, del_phi_2, cell_ptr_1, cell_ptr_2);

//     // cell_ptr_1->phi += del_phi_1;
//     // cell_ptr_2->phi += del_phi_2;
//     // del_F -= calculate_gradient_energy_dif_mesh(cell_ptr_2->phi, -del_phi_2, -del_phi_1, cell_ptr_2, cell_ptr_1);

//     // cell_ptr_1->phi -= del_phi_1;
//     // cell_ptr_2->phi -= del_phi_2;

//     return del_F;

// }

inline double Model::calculate_energy_dif(const Event* event_ptr){

    // IMPORTANT: Do not under any circumstance add in dependence on phi2 in here
    // Unless you adapt the selective updating process, the calculations would be wrong. 
    // This is because if an event targeting cell 2 is updated and not executed for many iterations,
    // the target concentration will likely have changed, so the propensity would be wrong.

    Cell *cell_ptr_1, *cell_ptr_2;
    double del_phi_1, del_phi_2;

    double del_f = 0;   // Total energy density
    double del_F = 0;   // Total energy

    cell_ptr_1 = &cells[event_ptr->target_cells[0]];
    cell_ptr_2 = &cells[event_ptr->target_cells[1]];

    del_phi_1 = -1/cell_ptr_1->capacity;   // Todo: change in next iteration to be particle dependent and cell size dependent
    del_phi_2 = 1/cell_ptr_2->capacity;

    // Enthalpy
    // cout << "Enthalpy difference cell "<< cell_ptr_1->idx << ": " << calculate_enthalpy_dif(cell_ptr_1->phi, del_phi_1) << endl;
    del_f += calculate_enthalpy_dif(cell_ptr_1->phi, del_phi_1);

    // cout << "Enthalpy difference cell "<< cell_ptr_1->idx << ": " << calculate_enthalpy_dif_adapted(cell_ptr_1->phi, cell_ptr_2->phi) << endl;
    // del_F += calculate_enthalpy_dif_adapted(cell_ptr_1->phi, cell_ptr_2->phi);

    // Surface tension
    // cout << "Surface difference cell "<< cell_ptr_1->idx << ": " << calculate_gradient_energy_dif_mesh(cell_ptr_1->phi, del_phi_1, del_phi_2, cell_ptr_1, cell_ptr_2) << endl;
    // del_f += calculate_surface_tension_dif_grid(cell_ptr_1->phi, del_phi_1, cell_ptr_1);
    del_F += calculate_gradient_energy_dif_mesh(cell_ptr_1->phi, del_phi_1, del_phi_2, cell_ptr_1, cell_ptr_2);
    
    // Repulsion
    // cout << "Repulsion difference cell "<< cell_ptr_1->idx << ": " << calculate_repulsion_dif(phi_1, del_phi_1) << endl;
    // del_f += calculate_repulsion_dif(cell_ptr_1->phi, del_phi_1);

    // cout << "Event " << event_ptr->idx << " in cell " << events[event_ptr->idx].target_cells[0] << endl;
    // cout << "Energy: " << del_f << endl;
    // cin.get();

    // if(fabs(del_f)>F_max){
    //     F_max = fabs(del_f);
    // }

    if(del_f*cell_ptr_1->size + del_F < -50 || false){
        cout << "Cell " << cell_ptr_1->idx << endl;
        cout << "Energy of individual step is " << del_f*cell_ptr_1->size + del_F << ", smaller than -50" << endl;
        cout << "Enthalpy difference cell "<< cell_ptr_1->idx << ": " << cell_ptr_1->size*calculate_enthalpy_dif(cell_ptr_1->phi, del_phi_1) << endl;
        cout << "Surface difference cell "<< cell_ptr_1->idx << ": " << calculate_gradient_energy_dif_mesh(cell_ptr_1->phi, del_phi_1, del_phi_2, cell_ptr_1, cell_ptr_2) << endl;
        // cout << "Repulsion difference cell "<< cell_ptr_1->idx << ": " << calculate_repulsion_dif(cell_ptr_1->phi, del_phi_1) << endl << endl;

        cout << "Density origin phi: " << cell_ptr_1->phi << ", new density phi + del_phi: " << cell_ptr_1->phi + del_phi_1 << endl;
        cin.get();

    }

    // If adding del_phi to cell 2 phi would cause going over capacity, energy is infinity;
    // if(phi_2 + del_phi_2 >= 1){
    //     del_f = numeric_limits<double>::infinity();
    // }

    // cout << "Gradient energy: " << del_F << endl;
    // cout << "Mixing energy: " << del_f*cell_ptr_1->size << endl;
    // cin.get();

    return del_f*cell_ptr_1->size + del_F;

}

inline double Model::calculate_enthalpy_dif(const double phi, const double del_phi){
    // Calculates enthalpy change for a single cell
    // Flory-Huggins mixing energy complete term
    // return -chi*del_phi*(2*phi - 1 + del_phi);

    // Flory-Huggins total energy if the only interactions are protein attractions
    // Going through a transition state with only solvent bound
    return -chi*del_phi*(2*phi + del_phi);
}

inline double Model::calculate_enthalpy_dif_adapted(const double phi_o, const double phi_t){
    // Calculates enthalpy change for a single cell
    // Flory-Huggins mixing energy complete term
    // return -chi*del_phi*(2*phi - 1 + del_phi);

    // Flory-Huggins total energy if the only interactions are protein attractions
    // Going through a transition state with only solvent bound
    if(phi_o > phi_t){
        return 2*chi*(phi_o-phi_t);
    } else{
        return 0;
    }
}


inline double Model::calculate_gradient_energy_dif_mesh(const double phi_o, const double del_phi_o, const double del_phi_t, Cell *cell_ptr_o, Cell *cell_ptr_t){
    // Calculates interfacial energy change for a single cell. See python file for correct formulas. Note that indices are confusing
    // Index 'o' is origin cell (i.e. cell 1), index 't' is target cell (i.e. cell 2 if in a shared face).
    // Index '2' and '3' refer to the other two cells. If a face is shared between origin and target, index 2 generally isn't the target cell, but may be.

    
    double del_F = 0;
    double del_F_per_face = 0;
    double phi_2, phi_3;

    unsigned int i_cell;

    for(Face *face : cell_ptr_o->adjacent_faces){
        
        // Energy times area for a given face. 
        del_F_per_face = 0;

        if(face->cells.size() != 3){
            cout << "Formula only works for triangles. If your faces are not triangles, god help you" << endl;
            cin.get(); exit(1);
        }

        // Look which cell is cell of change
        for(i_cell = 0; i_cell < 3; i_cell++){
            if (face->cells[i_cell] == cell_ptr_o){
                // Found i_cell that corresponds to cell 1
                break;
            }
        }

        // Densities of cells, we don't know if they are a target.
        phi_2 = face->cells[(i_cell+1)%3]->phi;
        phi_3 = face->cells[(i_cell+2)%3]->phi;


        del_F_per_face += 0.25*((2*phi_o*del_phi_o + del_phi_o*del_phi_o)*pow(face->lengths[i_cell], 2)
                        -2*face->theta_terms[(i_cell+1)%3]*phi_3*del_phi_o
                        -2*face->theta_terms[(i_cell+2)%3]*phi_2*del_phi_o);

        // Apply correction if face is shared with target. The correction is shared with the reverse reaction, so factor -2 becomes factor -1
        if (face->cells[(i_cell+1)%3] == cell_ptr_t){
            del_F_per_face += 0.25*(-face->theta_terms[(i_cell+2)%3]*del_phi_o*del_phi_t);
        } else if (face->cells[(i_cell+2)%3] == cell_ptr_t){
            del_F_per_face += 0.25*(-face->theta_terms[(i_cell+1)%3]*del_phi_o*del_phi_t);
        }

        // cout << "Face surface energy difference" << endl;
        // cout << "Origin cell: " << cell_ptr_o->idx << ", target cell: " << cell_ptr_t->idx << endl;
        // cout << "Cell 1, idx " << cell_ptr_o->idx << ": phi = " << phi_o << endl;
        // cout << "Cell 2, idx " << face->cells[(i_cell+1)%3]->idx << ": phi = " << phi_2 << endl;
        // cout << "Cell 3, idx " << face->cells[(i_cell+2)%3]->idx << ": phi = " << phi_3 << endl;
        // cout << "Energy difference for this face: " << del_F_per_face << endl;

        // cout << "Wrongly calculated energy: " << del_F_per_face << endl;
        // cout << "Correctly calculated energy: " << del_F_per_face/(face->area) << endl;
        // cout << "Factor difference: " << (face->area) << endl;
        // cin.get();

        del_F += del_F_per_face/(face->area);
        // del_F += del_F_per_face;
        
    }

    return del_F*kappa/2;

}

// inline double Model::calculate_surface_tension_dif_grid(const double phi, const double del_phi, Cell *cell_ptr){
//     // Calculates interfacial energy change for a single cell
//     Cell *neighbour_ptr;
//     double del_f = 0;
//     double phi_n;

//     for(unsigned int i_neighbour = 0; i_neighbour < cell_ptr->n_neighbours; i_neighbour++){
//         neighbour_ptr = &cells[cell_ptr->neighbours[i_neighbour]];

//         if(cell_ptr->dim == neighbour_ptr->dim){
//             phi_n = neighbour_ptr->phi;
//             del_f += (2*phi*del_phi + del_phi*del_phi -2*del_phi*phi_n) / pow(cell_ptr->neighbour_dist[neighbour_ptr->idx], 2);
//         }
//     }

//     return kappa*del_f;
// }

inline double Model::calculate_repulsion_dif(const double phi, const double del_phi){

    // Using Carnahan-Starling Free energy
    // return (4*pow(phi + del_phi, 2) - 3*pow(phi + del_phi, 3))/pow(1- phi - del_phi, 2) - (4*pow(phi, 2) - 3*pow(phi, 3))/pow(1- phi, 2);
    
    // Using infinitely high potential. Applied somewhere else.
    return 0;

}

inline void Model::update_diffusion_rate(Event* event_ptr, const double del_f){

    // Diffusion is from cell to neighbour
    Reaction *reaction_ptr = event_ptr->reaction_ptr;

    // Enthalpy factor only determines outward diffusion
    event_ptr->rate = exp(-del_f) * event_ptr->base_rate;
    // cout << "Diffusion rate without correction: " << (rate_constants[reaction_ptr->k_idx]) << endl;
    // cout << "Diffusion rate with distance correction: " << (rate_constants[reaction_ptr->k_idx]) / pow(cell_ptr->neighbour_dist[event_ptr->target_cells[1]], 2) << endl;
    // cout << "Diffusion rate with distance and energy correction: " << event_ptr->rate << endl;

    // cout << "Distance: " << cell_ptr->neighbour_dist[event_ptr->target_cells[1]] << endl;
    // cout << "Energy: " << del_f << endl << endl;

    if(event_ptr->rate < 0){
        cout << "Rate has gone negative" << endl;
        cin.get();
        exit(1);
    }

    if(event_ptr->rate > 1e50){
        cout << "Rate is very high (" << event_ptr->rate << ") because of energy " << del_f << endl;
        cout << "Event " << event_ptr->idx << ", Reaction " << reaction_ptr->reaction_name << endl;
        cin.get();
    }


}

inline double Model::get_event_propensity(const Event* event_ptr){

    Reaction *reaction_ptr;
    Cell *cell_ptr;

    double propensity;
    int reactant_number;
    int i_species, i_reaction_species;

    // // Quadratic terms
    // int state_change;

    reaction_ptr = event_ptr->reaction_ptr;

    // cout << "Target cells: " << endl;
    // for(int i_reactant = 0; i_reactant < reaction_ptr->n_reaction_species; i_reactant++){
    //     cout << "Cell " << event_ptr->target_cells[i_reactant] << endl;
    // }

    // cout << "Involved reactants and their number: " << endl;

    
    // if (event_ptr->idx == 14627 && i_iter > 700000){
    //     cout << "Iteration " << i_iter << endl;

    //     for(i_reaction_species = 0; i_reaction_species < reaction_ptr->n_reaction_species; i_reaction_species++){
    //         cell_ptr = &cells[event_ptr->target_cells[i_reaction_species]];

    //         cout << "New content for species: " << cell_ptr->content + reaction_ptr->state_change[i_reaction_species] << endl;
    //         cout << "Capacity: " << cell_ptr->capacity << endl;
    //         cout << "comparitson " << (cell_ptr->content + reaction_ptr->state_change[i_reaction_species] >= cell_ptr->capacity) << endl;
    //         cin.get();
    //     }

    // }


    // Check for a reaction that would cause capacity to overflow
    for(i_reaction_species = 0; i_reaction_species < reaction_ptr->n_reaction_species; i_reaction_species++){
        cell_ptr = &cells[event_ptr->target_cells[i_reaction_species]];

        if ((cell_ptr->content + reaction_ptr->state_change[i_reaction_species] >= cell_ptr->capacity)){

            if(reaction_ptr->state_change[i_reaction_species] < 0){
                cout << "Cell was already over capacity" << endl;
                cin.get();
            }
            // cout << "Rejected reaction based on capacity" << endl;
            // n_reject++;
            return 0;
        }
    }
    
    // TODO Optimize: double loop that kind of does the same thing
    
    propensity = event_ptr->rate;
    // cout << event_ptr->idx << " " << reaction_ptr->reaction_name << " in cell " << events[event_ptr->idx].target_cells[0] << endl;
    // cout << "Base rate: " << propensity << endl;

    for(int i_reactant = 0; i_reactant < reaction_ptr->n_reactants; i_reactant++){

        // Index "i_reactant" is an index over all reactants in the reaction.
        // It points to an entry in "reactants" which contains an index i_reaction_species.
        // This points to an entry in "species" in Model which is the index of the actual reaction.
        i_reaction_species = reaction_ptr->reactants[i_reactant];
        i_species = reaction_ptr->reaction_species[i_reaction_species];
        
        
        // Look at target cell, and in that cell, get number of molecules of the species.
        reactant_number = cells[event_ptr->target_cells[i_reaction_species]].state[i_species];
        propensity *= double(reactant_number);

        // Quadratic terms. Commented out to speed up the model.
        // state_change = reaction_ptr->state_change[i_reaction_species];

        // For quadratic terms, the c values are twice as high as the k values. (See Gillespie, 1977).
        // Moreover, there are 1/2X(X-1) pairs instead of XY. These effects cancel, so we obtain:
        // if(state_change == -2){
        //     propensity *= (reactant_number-1);
        // }
        // else if(state_change != -1 || state_change != -2){
        //     cout << "State change " << state_change << "is not linear or quadratic, and hence not covered by model" << endl;
        //     cin.get();
        //     exit(1);
        // }

        // cout << i_species << " " << species[i_species] << ": ";
        // cout << reactant_number << endl;
        // cin.get();

    }

    if(propensity < 0){
        cout << "Propensity for event " << event_ptr->idx << " went negative" << endl;
        cin.get(); exit(1);
    }

    if(propensity > 1e50){
        cout << "Propensity for event " << event_ptr->idx << " is very large" << endl;
        cin.get(); exit(1);
    }

    return propensity;

}

inline void Model::recalculate_propensities(){

    int i_cell;
    double propensity;
    double del_f;

    set<int> visited_cells;
    Reaction *reaction_ptr_orig;
    Event *event_ptr, *event_ptr_orig;
    Cell *cell_ptr;
    // Cell *neighbour_ptr;

    // TODO: make function update_event that is called. Make it return propensity, not update it directly

    event_ptr_orig = &events[mu];
    reaction_ptr_orig = event_ptr_orig->reaction_ptr;

    // cout << reaction_ptr_orig->reaction_name << " is executed" << endl;
    // cout << "Number of cells to update: " << reaction_ptr_orig->n_reaction_species << endl;

    // cout << "Target cells: " << endl;
    // for(int i_reactant = 0; i_reactant < reaction_ptr_orig->n_reaction_species; i_reactant++){
    //     cout << "Cell " << event_ptr_orig->target_cells[i_reactant] << endl;
    // }

    for(int i_target_cell = 0; i_target_cell < reaction_ptr_orig->n_reaction_species; i_target_cell++){

        i_cell = event_ptr_orig->target_cells[i_target_cell];

        // if not visited this cell already, visit it
        if(visited_cells.find(i_cell) == visited_cells.end()){

            // cout << "Entered in target cell loop for cell affected by reaction species " << i_target_cell << " in cell " << i_cell << endl;

            cell_ptr = &cells[i_cell];

            t1 = timing_clock.now();

            // Loop over all affected events, and update them. This includes only all directly affected cells
            for(unsigned int i_affected_event = 0; i_affected_event < cell_ptr->n_affected_events; i_affected_event++){
                event_ptr = cell_ptr->affected_events[i_affected_event];

                // cout << "Updating affected event " << event_ptr->reaction_ptr->reaction_name << " in cell(s)";
                // for(int i_reaction_species = 0; i_reaction_species < event_ptr->reaction_ptr->n_reaction_species; i_reaction_species++){
                //     cout << " " << event_ptr->target_cells[i_reaction_species];
                // } cout << endl;

                // Update diffusion events
                if(event_ptr->reaction_ptr->is_diffusion_like && phase_separate) {
                    // cout << "Updating phase seperation parameters" << endl;
                    del_f = calculate_energy_dif(event_ptr);
                    update_diffusion_rate(event_ptr, del_f);
                }

                // Update the propensity at that location
                propensity = get_event_propensity(event_ptr);

                // cout << "Original propensity of event " << event_ptr->idx << ": " << propensities[event_ptr->idx] << endl;
                // cout << "Original total propensity " << a0 << endl;

                a0 -= propensities[event_ptr->idx];
                propensities[event_ptr->idx] = propensity;
                a0 += propensity;

                // cout << "New propensity of event " << event_ptr->idx << ": " << propensities[event_ptr->idx] << endl;
                // cout << "New total propensity " << a0 << endl;
                // cin.get();
            }

            // cin.get();

            t2 = timing_clock.now();
            affected_cell_time += duration_cast<milliseconds>(t2 - t1);


            // Old comment, to be moved.
            // cout << "Going into neighbour checking loop" << endl;
            // Checking that any neighbouring events still do not cause capacity overflow. E.g. can be the case if previously, a cell was not near capacity,
            // so propensity for diffusion into this cell is nonzero.
            // After updates, cell can be at capacity, but propensities of other events would not have changed. This loop is therefore needed.

            visited_cells.insert(i_cell);
        }
    }

}

inline void Model::calculate_propensities(){

    Event* event_ptr;
    Cell* cell_ptr;

    double del_f;
    double propensity;

    // Reset total propensity
    a0 = 0;

    // Loop over all cells to determine diffusion correction factor
    for(int i_cell = 0; i_cell < n_cells; i_cell++){
        cell_ptr = &cells[i_cell];

        cell_ptr->content = 0;
        cell_ptr->phi = 0;

        // Count all particles in cells
        for(int i_species = 0; i_species < n_species; i_species++){
            cell_ptr->content += cell_ptr->state[i_species];
            cell_ptr->phi = cell_ptr->content/cell_ptr->capacity;
        }
    }


    // Loop over all events, and calculate propensities for each event.
    for(int i_event = 0; i_event < n_events; i_event++){

        event_ptr = &events[i_event];
        
        // Update diffusion events
        if(event_ptr->reaction_ptr->is_diffusion_like && phase_separate) {
            del_f = calculate_energy_dif(event_ptr);
            update_diffusion_rate(event_ptr, del_f);
        }

        propensity = get_event_propensity(event_ptr);

        propensities[i_event] = propensity;
        a0 += propensity;   

        // if(verbose){
        //     cout << "Propensity: " << propensity << endl;
        //     cin.get();
        // }
    }
}

inline void Model::update_state(){

    Reaction *reaction_ptr;
    Event *event_ptr;
    Cell *cell_ptr;

    event_ptr = &events[mu];
    reaction_ptr = event_ptr->reaction_ptr;

    // Debugging code
    // cout << event_ptr->idx << " " << reaction_ptr->reaction_name << " in cell " << event_ptr->target_cells[0];
    // if (reaction_ptr->is_diffusion_like){
    //     cout << " to cell " << event_ptr->target_cells[1];
    // }
    // cout << ", a = " << propensities[mu] << "   " << endl;


    // Count the reaction
    reaction_ptr->count++;

    for(int i_reaction_species = 0; i_reaction_species < reaction_ptr->n_reaction_species; i_reaction_species++){
        
        // Get pointer to target cell
        cell_ptr = &(cells[event_ptr->target_cells[i_reaction_species]]);

        // Update state of target cell
        cell_ptr->state[reaction_ptr->reaction_species[i_reaction_species]] += reaction_ptr->state_change[i_reaction_species];
        
        // Update number of particles in the cell
        cell_ptr->content += reaction_ptr->state_change[i_reaction_species];
        cell_ptr->phi = cell_ptr->content/cell_ptr->capacity;

        // Safety check: population of state must always be 0 or bigger
        if(cell_ptr->state[reaction_ptr->reaction_species[i_reaction_species]] < 0){
            cout << "Population of cell  " << cell_ptr->idx << ", reactant " << i_reaction_species << " went negative." << endl;
            cout << "Last reaction: " << events[mu].reaction_ptr->reaction_name << endl;
            error_dump();
            cin.get(); exit(1);
        }

        if(cell_ptr->content > cell_ptr->capacity){
            cout << "Content of cell " << cell_ptr->idx <<  " exceeded capacity." << endl;
            cout << "Last reaction: " << events[mu].reaction_ptr->reaction_name << endl;

            cout << cell_ptr->content << " " << cell_ptr->capacity << endl;
            cout << "Event " << event_ptr->idx << endl;

            error_dump();
            cin.get(); exit(1);
        }
    }
}

void Model::add_rate_constant_sweep(const int i_rate_constant, const double k_start, const double k_end, const int n_points, const string sweep_type){

    if(n_points <= 1){
        tell("Number of points in the sweep must be bigger than 1", 1);
    }

    vector<double> sweep;

    if(sweep_type == "lin"){
        for(int i_point = 0; i_point < n_points; i_point++){
            sweep.push_back(i_point*(k_end-k_start)/(n_points-1) + k_start);
        }
    }
    else if(sweep_type == "log"){

        if(k_start == 0){
            tell("In log sweep, start number cannot be 0", 1);
        }

        for(int i_point = 0; i_point < n_points; i_point++){
            sweep.push_back(k_start * pow(double(k_end)/k_start, double(i_point)/(n_points-1)));
        }
    }

    // Add to model member
    sweeps_rate_constants.push_back(sweep);
    sweeps_rate_constant_idx.push_back(i_rate_constant);

}


double distance(vector<double> v1, vector<double> v2){

    double dist2 = 0;

    for(unsigned int i = 0; i < v1.size(); i++){
        dist2 += pow(v1[i] - v2[i], 2);
    }

    return sqrt(dist2);

}

void Model::create_spherical_grid(const int n_points, const double r){

    if(n_points > MAX_CELLS){
        tell("Number of cells exceeds maximum number of cells", 1);
    }

    
    string filename = "./Grids/" + to_string(n_points) + ".txt";

    ifstream point_file;
    point_file.open(filename);
    if (!point_file){
        cout << "File " << filename << " does not exist" << endl;
        cin.get();
        exit(1);
    }

    n_cells = n_points;
    n_cells_per_dim[2] = n_cells;
    n_cells_per_dim[3] = 1;

    // Bulk component
    cells[n_cells].idx = n_cells;
    cells[n_cells].pos.x = 0;
    cells[n_cells].pos.y = 0;
    cells[n_cells].pos.z = 0;

    cells[n_cells].size     = 4./3.*M_PI*pow(r, 3);
    cells[n_cells].capacity = cells[n_cells].size * particles_per_volume;
    cells[n_cells].dim = 3;

    cells[n_cells].neighbours.resize(n_cells);
    cells[n_cells].affected_events.resize(n_cells*MAX_REACTIONS);

    string line;
    vector<vector<double>> coordinates(n_points, vector<double>(3, 0));

    for(int i_cell = 0; i_cell < n_cells; i_cell++) {

        for(int i = 0; i < 3; i++){
            getline(point_file, line);
            coordinates[i_cell][i] = stod(line);
        }

        cells[i_cell].pos.x = r*coordinates[i_cell][0];
        cells[i_cell].pos.y = r*coordinates[i_cell][1];
        cells[i_cell].pos.z = r*coordinates[i_cell][2];

        cells[i_cell].idx = i_cell;
        cells[i_cell].dim = 2;

    }

    // for (auto & x : coordinates) {
    //     cout << x[0] << x[1] << x[2] << endl;
    // }

    double neighbour_dist;
    double min_dist = cells[0].distance(cells[1]);

    // Finding minimal distance
    for(int cell_idx = 0; cell_idx < n_cells; cell_idx++){
        for(int neighbour_idx = cell_idx + 1; neighbour_idx < n_cells; neighbour_idx++){

            neighbour_dist = cells[cell_idx].distance(cells[neighbour_idx]);
            if (neighbour_dist < min_dist){min_dist = neighbour_dist;}
            // distances.push_back(distance(coordinates[cell_idx], coordinates[neighbour_idx]));

        }
    }
    cout << "Creating spherical grid" << endl;
    cout << "Minimal distance: " << min_dist << endl;

    vector<Edge> edges;
    int n_edges;

    // Adding neighbours
    for(int cell_idx = 0; cell_idx < n_cells; cell_idx++){
        for(int neighbour_idx = cell_idx + 1; neighbour_idx < n_cells; neighbour_idx++){
            if (cell_idx == neighbour_idx) {continue;}
            neighbour_dist = cells[cell_idx].distance(cells[neighbour_idx]);
            if(neighbour_dist < 1.35*min_dist){
                cells[cell_idx].add_neighbour(neighbour_idx, neighbour_dist, 0);
                cells[neighbour_idx].add_neighbour(cell_idx, neighbour_dist, 0);

                // Adding edge for face finding
                Edge edge(cells[cell_idx], cells[neighbour_idx]);
                edges.push_back(edge);
                // The interface is set at 0. This is because when deciding neighbours, the interface is not yet known
                // We calculate the interface later, and set it.
                // cout << cells[cell_idx].distance(cells[neighbour_idx]) << endl;
            }
        }
    }

    n_edges = edges.size();

    // Finding faces
    // for i in range(0, len(edges) - 1):
    // for j in range(i + 1, len(edges)):
    //     e1 = edges[i]
    //     e2 = edges[j]
    //     if e1[0] == e2[0]:
    //         if (e1[1], e2[1]) in edges:
    //             print(e1[0] + e1[1] + e2[1])
    //     else:
    //         break


    // Finding the faces. Each face is a triangle, and has a list of cells that form the vertices of the triangle.
    // Subsequently, some geometric quantities of the triangle area calculated.

    // If you don't reserve, weird shit happens. References go awry
    // Use Eulers characteristic: V - E + F = 2. Rewrite to F = E - V + 2
    faces.reserve(n_edges - n_points + 2);

    for(unsigned int i_edge = 0; i_edge < edges.size(); i_edge++){
        for(unsigned int j_edge = i_edge + 1; j_edge < edges.size(); j_edge++){

            // Finding common first vertex
            if (edges[i_edge].c1 == edges[j_edge].c1){
                
                // Finding bridging edge that complete triangle
                for (unsigned int k_edge = 0; k_edge < edges.size(); k_edge++){
                    if ((edges[i_edge].c2 == edges[k_edge].c1) && (edges[j_edge].c2 == edges[k_edge].c2)){
                        
   
                        // Generate face and save in model
                        Face face(edges[i_edge].c1, edges[i_edge].c2, edges[j_edge].c2);
                        face.generate_angles_and_lengths();
                        faces.push_back(face);
                        // cout << "Triangle " << face.cells[0]->idx << ", " << face.cells[1]->idx << ", " << face.cells[2]->idx << endl;

                        // Save reference in cell
                        for(Cell *cell : face.cells){
                            cell->adjacent_faces.push_back(&faces.back());
                        }
                    }
                }
            } else {
                break;
            }
        }
    }

    // Check that the area of the faces are calculated correctly.
    // double tot_area = 0;
    // for (Face face : faces){
    //     tot_area += face.area;

    // }
    // cout << "Total area: " << tot_area << endl;
    // cin.get();

    // for (Face *face : cells[512].adjacent_faces){
    //     cout << "Triangle " << face->cells[0]->idx << ", " << face->cells[1]->idx << ", " << face->cells[2]->idx << endl;
    //     cout << "First normal " << face->edge_normals[0].abs() << endl;
    // }

    // Finding loops around a centre cell
    int start_cell_idx, centre_cell_idx;
    Cell *centre_cell_ptr, *cur_cell_ptr, *neighbour_ptr;
    set<int> visited_cells;

    for(int cell_idx = 0; cell_idx < n_cells; cell_idx++){

        int n_loop = 1;
        visited_cells.clear();

        // Finding centre cell and saving index
        centre_cell_ptr = &cells[cell_idx];
        centre_cell_idx = centre_cell_ptr->idx;

        // Starting circle around centre at random neighbour cell, and saving this cell.
        cur_cell_ptr = &cells[centre_cell_ptr->neighbours[0]];
        start_cell_idx = cur_cell_ptr->idx;

        // Save which neighbour the vertex belongs too
        
        centre_cell_ptr->ordered_neighbours.push_back(cur_cell_ptr->idx);
        visited_cells.insert(cur_cell_ptr->idx);

        // cout << "Starting from centre " << centre_cell_idx << endl;
        // cout << "Starting loop" << endl;
        // cout << cur_cell_ptr->idx << endl;

        // While loop that stops when the current cells (again) neigbours the starting cell.
        do {
            for(unsigned int i_neighbour = 0; i_neighbour < cur_cell_ptr->n_neighbours; i_neighbour++){
                neighbour_ptr = &cells[cur_cell_ptr->neighbours[i_neighbour]];
                // If potential new cell in chain also borders the centre cell, and has not yet been visited, go to this cell

                // cout << "In loop for neighbour cell " << neighbour_ptr->idx << endl;
                // cout << neighbour_ptr->has_neighbour(centre_cell_idx) << endl;
                // cout << (visited_cells.find(neighbour_ptr->idx) == visited_cells.end()) << endl;

                

                if(neighbour_ptr->has_neighbour(centre_cell_idx) && (visited_cells.find(neighbour_ptr->idx) == visited_cells.end())){

                    // cout << "Neighbour cell " << neighbour_ptr->idx << " satisfies" << endl;

                    // Find centroid in here, and save to centre cell
                    Vertex centroid = (centre_cell_ptr->pos + cur_cell_ptr->pos + neighbour_ptr->pos)/3;
                    centre_cell_ptr->vertices.push_back(centroid);

                    // Update current cell to find next cell
                    cur_cell_ptr = neighbour_ptr;
                    visited_cells.insert(cur_cell_ptr->idx);
                    
                    // Save which neighbour the vertex belongs too
                    centre_cell_ptr->ordered_neighbours.push_back(cur_cell_ptr->idx);

                    n_loop++;
                    break;
                }
            }

            if(n_loop == 1){
                cout << "No neighbour also bordering centre cell " << centre_cell_idx << " found. Loop may not exist" << endl;
                cin.get();
                exit(1);
            }

        } while(!cur_cell_ptr->has_neighbour(start_cell_idx) || n_loop == 2);
        
        // Adding first ordered neighbour also as last ordered neighbour
        centre_cell_ptr->ordered_neighbours.push_back(centre_cell_ptr->ordered_neighbours.front());

        // Adding centroid between final and starting cell
        Vertex centroid = (centre_cell_ptr->pos + cur_cell_ptr->pos + cells[start_cell_idx].pos)/3;
        centre_cell_ptr->vertices.push_back(centroid);

        // Adding last centroid at the begining
        centre_cell_ptr->vertices.insert(centre_cell_ptr->vertices.begin(), centre_cell_ptr->vertices.back());

        // Setting the width of the interfaces between two cells of dual grid
        for (unsigned int i_neighbour = 0; i_neighbour < centre_cell_ptr->n_neighbours; i_neighbour++){
            // cout << "Distance:" << centre_cell_ptr->vertices[i_neighbour].distance(centre_cell_ptr->vertices[i_neighbour + 1]) << endl;
            // 1. Setting interface width between cells
            centre_cell_ptr->neighbour_interface[centre_cell_ptr->ordered_neighbours[i_neighbour]] = centre_cell_ptr->vertices[i_neighbour].distance(centre_cell_ptr->vertices[i_neighbour + 1]);
        }

        // Calculating areas via cross products or something similar
        Vertex area_vector(0, 0, 0);
        for (unsigned int i_neighbour = 1; i_neighbour < centre_cell_ptr->n_neighbours - 1; i_neighbour++){
            area_vector = area_vector + (centre_cell_ptr->vertices[i_neighbour] - centre_cell_ptr->vertices[0]).cross(
                                         centre_cell_ptr->vertices[i_neighbour + 1] - centre_cell_ptr->vertices[0]);
        }
        
        centre_cell_ptr->size = 0.5*area_vector.abs();
        centre_cell_ptr->capacity = centre_cell_ptr->size * particles_per_area;
        // cout << "Area of cell " << centre_cell_ptr->idx << ": " << centre_cell_ptr->size << endl;
        // cout << "Particles per area: " << particles_per_area << endl;
        // cin.get();

        // Adding bulk component neighbour
        cells[cell_idx].add_neighbour_alt_dim(n_cells, cells[cell_idx].size);
        cells[n_cells].add_neighbour_alt_dim(cell_idx, cells[cell_idx].size);

        // for (Vertex vertex : centre_cell_ptr->vertices){
        //     cout << "Vertex " << vertex.x << " " << vertex.y << " " << vertex.z << endl;
        //     cout << "Length " << vertex.abs() << endl;
        // }

        // for (int idx : centre_cell_ptr->ordered_neighbours){
        //     cout << "Neighbour " << idx << endl;
        // }

        // cin.get();
    }

    // Update number of cells for bulk component
    n_cells++;

    // // Check that the area is calculated correctly
    // double tot_area = 0;
    // for(int cell_idx = 0; cell_idx < n_cells; cell_idx++){
    //     if (cells[cell_idx].dim == 2){
    //         tot_area += cells[cell_idx].size;
    //     }
    // }
    // cout << "Total area: " << tot_area << endl;
    // cin.get();

    cout << "Number of cells: " << n_cells << endl;
}


void Model::create_cartesian_grid(const double d_x, const double d_y, const double d_z, const int n_x, const int n_y, const bool periodic){

    n_cells = 0;

    if(n_x*n_y > MAX_CELLS){
        tell("Number of cells exceeds maximum number of cells", 1);
    }
    if(n_x <= 1 || n_y <= 1){
        tell("Grid must have a 2 or more points in both dimensions", 1);
    }

    if(d_x <= 0 || d_y <= 0 || d_z <= 0){
        tell("All lengths must be nonzero and positive", 1);
    }

    double dist_x, dist_y;
    dist_x = d_x/(n_x-1);
    dist_y = d_y/(n_y-1);

    n_cells = n_x*n_y;
    n_cells_per_dim[2] = n_cells;
    n_cells_per_dim[3] = 1;

    // Bulk component
    cells[n_cells].idx = n_cells;
    cells[n_cells].pos.x = 0;
    cells[n_cells].pos.y = 0;
    cells[n_cells].pos.z = 0;

    cells[n_cells].size     = d_x*d_y*d_z;
    cells[n_cells].capacity = cells[n_cells].size * particles_per_volume;
    cells[n_cells].dim = 3;

    cells[n_cells].neighbours.resize(n_cells);
    cells[n_cells].affected_events.resize(n_cells*MAX_REACTIONS);
    cout << "Cell capacity bulk component: " << cells[n_cells].capacity << endl;
    cout << "Cell capacity membrane: " << (d_x/n_x)*(d_y/n_y)*particles_per_volume << endl;

    // Creating periodic boundary conditions
    for(int i = 0; i < n_cells; i++){
        cells[i].idx = i;
        cells[i].pos.x = double(i%n_x)*dist_x;
        cells[i].pos.y = double(i/n_y)*dist_y;
        cells[i].pos.z = 0;

        cells[i].size   = (d_x/n_x)*(d_y/n_y);
        cells[i].capacity = cells[i].size * particles_per_area;

        cells[i].dim = 2;

        if(i%n_x != 0) {        // Cells not at left boundary
			cells[i].add_neighbour(i - 1, dist_x, d_y/n_y);
		} else if(periodic) {   // Cells at left boundary
			cells[i].add_neighbour(i + n_x - 1, dist_x, d_y/n_y);
		}
		
		if(i%n_x != n_x - 1) {  // Cells not at right boundary
			cells[i].add_neighbour(i + 1, dist_x, d_y/n_y);
		} else if(periodic) {   // Cells at right boundary
			cells[i].add_neighbour(i - n_x + 1, dist_x, d_y/n_y);
		}
		
		if(i > n_x - 1) {       // Cells not at top boundary
			cells[i].add_neighbour(i - n_x, dist_y, d_x/n_x);
		} else if(periodic) {   // Cells at top boundary
			cells[i].add_neighbour(i + n_cells - n_x, dist_y, d_x/n_x);
		}
		
		if(i < n_cells - n_x) { // Cells not at bottom boundary
			cells[i].add_neighbour(i + n_x, dist_y, d_x/n_x);
		} else if(periodic) {   // Cells at bottom boundary
			cells[i].add_neighbour(i%n_x, dist_y, d_x/n_x);
		}

        // Adding neighbour to bulk component
        cells[i].add_neighbour_alt_dim(n_cells, cells[i].size);
        cells[n_cells].add_neighbour_alt_dim(i, cells[i].size);

    }

    // Update number of cells for bulk component
    n_cells++;
}

void Model::load_phase_separation_parameters(){


    if(phase_separate){
        if(rate_constant_map.count("chi")){
            chi = rate_constants[rate_constant_map["chi"]];
            cout << "Assigned internal parameter chi = " << chi << endl;
        } else{
            cout << "Parameter chi is not defined" << endl;
            cin.get(); exit(1);
        }

        if(rate_constant_map.count("kappa")){
            kappa = rate_constants[rate_constant_map["kappa"]];
            cout << "Assigned internal parameter kappa = " << kappa << endl;
        } else{
            cout << "Parameter kappa is not defined" << endl;
            cin.get(); exit(1);
        }

        if(rate_constant_map.count("particles_per_area")){
            particles_per_area = rate_constants[rate_constant_map["particles_per_area"]];
            cout << "Assigned internal parameter particle per area = " << particles_per_area << endl;
        } else{
            cout << "Amount of particles per area is not defined" << endl;
            cin.get(); exit(1);
        }

        if(rate_constant_map.count("particles_per_volume")){
            particles_per_volume = rate_constants[rate_constant_map["particles_per_volume"]];
            cout << "Assigned internal parameter particle per volume = " << particles_per_volume << endl;
        } else{
            cout << "Amount of particles per volume not defined" << endl;
            cin.get(); exit(1);
        }

    }

    if(crowding){
        if(rate_constant_map.count("free_volume_fraction")){
            free_volume_fraction = rate_constants[rate_constant_map["free_volume_fraction"]];
            free_size_fraction[3] = free_volume_fraction;
        } else{
            cout << "Free volume fraction is not defined" << endl;
            cin.get(); exit(1);
        }

        if(rate_constant_map.count("free_area_fraction")){
            free_area_fraction = rate_constants[rate_constant_map["free_area_fraction"]];
            free_size_fraction[2] = free_area_fraction;
        } else{
            cout << "Free area fraction is not defined" << endl;
            cin.get(); exit(1);
        }
    }



}


void Model::run(){

    t_start = timing_clock.now();
    
    // Define events
    // For all unique reactions (actual reactions and diffusion), events are generated
    // so that for each cell, diffusions in all directions and reactions are stored in events.

    // Define chi and kappa based on input files.
    load_phase_separation_parameters();
    
    // Reset n_events in case that we do parameter sweeps. We have to redefine our events.
    reset();

    // Check if number of cells is small, affects validity of some calculations
    // since I made the assumption that one cell has unique neighbours. In a periodix 2x2 grid,
    // each cell only has 2 unique neighbours, not 4.
    if(n_cells < 10){
        cout << "Number of cells is small. If one cell can border a neighbouring cell more than once (e.g. in a periodic 2x2 grid) the energy difference formula fails" << endl;
        cin.get();
        cout << "Continue at own risk" << endl;
        cin.get();
    }

    // TODO: TODO: TODO: add check so that cells will NEVER have more contents than capacity

    // Diffusion reactions: diffusion of one species between neighbouring cells of the same dimension
    define_diffusion_reactions();

    // Normal reactions: reactions of multiple species in a single cell
    define_normal_reactions();

    // Boundary reactions: reactions of multiple species in multiple cells of different dimensions
    // E.g. attachment of a species from cytosol to the membrane.
    define_boundary_reactions();

    cout << ". . . . . . . . . . ."<< endl;
    cout << "Number of normal reactions: " << n_reactions << endl;
    cout << "Number of boundary reactions: " << n_boundary_reactions << endl;
    cout << "Number of species: " << n_species << endl;
    cout << "Number of cells: " << n_cells << endl;
    cout << "Number of events: " << n_events << endl;
    cout << "Done setting up events. Starting run" << endl;
    cout << ". . . . . . . . . . ."<< endl;


    // Recording initial conditions
    record_parameters();
    record_initial_conditions();

    // Recording number of cells and their neighbours
    record_cell_data();

    // Recording event data
    record_event_data();

    double a1 = 0, a_stable;
    double r_t, dt;

    Reaction *reaction_ptr;

    cout << ". . . . . . . . . . ."<< endl;
    cout << "Starting main loop" << endl;
    cout << "Kappa: " << kappa << endl;
    cout << "Chi: " << chi << endl;
    cout << ". . . . . . . . . . ."<< endl;

    // Main loop
    //while((i_iter < n_iterations) && t < t_max){
    while((i_iter < n_iterations)){
        // Recording concentrations
        if(i_iter%break_iter == 0){	
            cout << ". . . . . . . . . . ."<< endl;
            cout << "Iteration " << i_iter << endl;
            //cout << "Max energy: " << F_max << endl;
            //cout << "Number of rejected reactions: " << n_reject << endl;
            //cout << "Number of recalculations: " << n_recalculate << endl;
            //cout << "Previous stable propensity: " << a_stable << endl;

            // Find current time and convert to string
            t_now = timing_clock.now();
            total_time = duration_cast<seconds>(t_now - t_start);
            //cout << "Seconds since start of run: " << total_time.count() << endl << endl;

            // Instability check. If total propensity no longer matches sum of individual propensities
            // due to numerical instability, wrong events can be selected.
            a1 = 0;
            for (int i_event = 0; i_event < n_events; i_event++) {
                a1 += propensities[i_event];
            }

            if (fabs(fabs(a0/a1) - 1) > 1e-8){
                cout << "Iteration " << i_iter << endl;
                cout << "Instability factor: " << fabs(fabs(a0/a1) - 1) << endl; 
                cout << "a0:" << a0 << endl;
                cout << "a1:" << a1 << endl;
                // cin.get();
            }
            // End instability check

            // Resets floating point error which may have accumulated using recalculate propensities
            // Save total propensity, if total propensity changes magnitude too much, calculate everything again
            calculate_propensities();
            a_stable = a0;

            t1 = timing_clock.now();
            record_concentrations(i_iter);
            t2 = timing_clock.now();
            record_time += duration_cast<milliseconds>(t2 - t1);
        }


        // Make sure that numerical accuracy doesn't diverge too quickly.
        if(fabs(log(a_stable) - log(a0)) > 5){
            calculate_propensities();
            a_stable = a0;
            n_recalculate++;
        }

        // // Cdc42 test.
        // const_count = 0;
        // for(int i_cell = 0; i_cell < n_cells; i_cell++){
        //     const_count += cells[i_cell].state[0];
        //     const_count += cells[i_cell].state[3];
        //     const_count += cells[i_cell].state[4];
        //     const_count += cells[i_cell].state[6];

        // }
        // cout << "Total cdc42: " << const_count << endl;

        // if (const_count != 1000){
        //     cout << "Event: " << events[mu].reaction_ptr->reaction_name << endl;
        //     cin.get();
        // }

        // // Bem1 test.
        // const_count = 0;
        // for(int i_cell = 0; i_cell < n_cells; i_cell++){
        //     const_count += cells[i_cell].state[1];
        //     const_count += cells[i_cell].state[7];
        //     const_count += cells[i_cell].state[8];

        // }
        // cout << "Total Bem1: " << const_count << endl;

        // GEF test.
        // const_count = 0;
        // for(int i_cell = 0; i_cell < n_cells; i_cell++){
        //     const_count += cells[i_cell].state[2];
        //     const_count += cells[i_cell].state[8];
        //     const_count += cells[i_cell].state[9];

        // }
        // cout << "Total GEF: " << const_count << endl;


        // Picking time of reaction	
        // Evaluates the time before the next reaction will take place, storing it in dt (Eq. 21a in Gillespie 1977)
        
        // Starting time measurement
        t1 = timing_clock.now();

        // Safety check on a0 not getting too small
        if(a0 < 1e-8){
            cout << "Total propensity has become too small" << endl;
            error_dump();
            cin.get();
            exit(1);
        }

        r_t = dis(gen);
        dt = -log(r_t)/a0;
        t += dt;
        
        // Picking event, i.e. picks the reaction that is going to take place (Eq. 21b in Gillespie 1977)
        a1 = 0;
        r_mu = dis(gen);

        for (mu = 0; mu < n_events; mu++) {

            a1 += propensities[mu];

            if (r_mu*a0 <= a1){
                break;
            }
        }
        
        
        // ---------------------------------------------------------------------------------
        //                                 Safety checks
        // ---------------------------------------------------------------------------------
        // Check that valid event is selected, i.e. not that the previous loop did not break.
        if (mu == n_events){
            cout << "No event was selected at iteration " << i_iter << ". The sum of all propensities no longer matches total propensity." << endl;
            cout << "This likely occured due to numerical instability: if propensities grow too large or if propensities are infinity." << endl;
            cin.get();
            exit(1);
        }

        // Save mu for logging
        mu_history.push_back(mu);
        
        // Ending time measurement and calculating time
        t2 = timing_clock.now();
        find_mu_time += duration_cast<milliseconds>(t2 - t1);

        // Timing and other stats
        if((i_iter%break_iter == 0) && timing){
            //cout << "Timing measures" << endl;
            //cout << "Time to choose mu: " << find_mu_time.count() << endl;
            //cout << "Time for directly affected cells: " << affected_cell_time.count() << endl;
            //cout << "Time for neighbour cells: " << neighbour_cell_time.count() << endl;
            //cout << "Time to update state: " << update_time.count() << endl;
            //cout << "Time to record: " << record_time.count() << endl << endl;;

            // // Most recent executed event
            // cout << "Most recent event: " << events[mu].idx << endl;
            // cout << events[mu].reaction_ptr->reaction_name << " in cell " << events[mu].target_cells[0];
            // cout << " (phi " << cells[events[mu].target_cells[0]].content << ")";
            // if (events[mu].reaction_ptr->is_diffusion_like){
            //     cout << " to cell " << events[mu].target_cells[1];
            //     cout << " (phi " << cells[events[mu].target_cells[1]].content << ")";
            // }
            // cout << ", a = " << propensities[mu] << "   " << endl;

        }

        // Verbose statements for debugging
        if((i_iter%break_iter == 0) && verbose){

            for(int i_event = 0; i_event < n_events; i_event++){
                
                reaction_ptr = events[i_event].reaction_ptr;
                cout << i_event << " " << reaction_ptr->reaction_name << " in cell(s)";
                for(int i_reaction_species = 0; i_reaction_species < reaction_ptr->n_reaction_species; i_reaction_species++){
                    cout << " " << events[i_event].target_cells[i_reaction_species];
                }
                cout << ", a = " << propensities[i_event] << "   " << endl;
            }

            cout << "a0 = "<< a0 << endl;
            cout << "Mu = "<< mu << endl;

            for(int i_cell = 0; i_cell < n_cells; i_cell++){
                for(int i_species = 0; i_species < n_species; i_species++){
                    cout << cells[i_cell].state[i_species] << "   ";
                }
                cout << endl;
            }
            cin.get();
            
        }


        // Updating the state vector with the chosen reaction
        t1 = timing_clock.now();
        update_state();
        t2 = timing_clock.now();
        update_time += duration_cast<milliseconds>(t2 - t1);

        // Recalculate propensities; only update cells affected by reaction mu.
        recalculate_propensities();

        // Update iteration number
        i_iter++;

    }

    record_concentrations(i_iter);
    cout << "Run " << i_run << " finished" << endl << endl << endl << endl;
}


void Model::do_parameter_sweep(){

    int n_particles, i_species;

    int i_rate_constant;
    double rate;

    // Internal variable for looping over all combinations of parameters.
    int i_parameter;    

    // Count how many different configurations are possible

    // Repeats
    int n_model_configs = n_repeats;

    // Species
    for(unsigned int i_vector = 0; i_vector < sweeps_species_number.size(); i_vector++){
        n_model_configs *= sweeps_species_number[i_vector].size();
    }
    // Rate constants
    for(unsigned int i_vector = 0; i_vector < sweeps_rate_constants.size(); i_vector++){
        n_model_configs *= sweeps_rate_constants[i_vector].size();
    }

    for(i_run = 0; i_run < n_model_configs; i_run++){
        i_parameter = i_run/n_repeats;
        cout << "Run " << i_run << endl;
        cout << "Parameter combination " << i_parameter << endl;
        cout << "Repeat " << i_run%n_repeats << endl;


        cout << "Size of sweeps:" << endl;
        cout << "Rate constants:" << sweeps_rate_constants.size() << endl;
        cout << "Species:" << sweeps_species_number.size() << endl;

        // Rate constants
        for(unsigned int i_vector = 0; i_vector < sweeps_rate_constants.size(); i_vector++){

            // Find number of particles and distribute them
            rate = sweeps_rate_constants[i_vector][i_parameter%sweeps_rate_constants[i_vector].size()];
            i_rate_constant = sweeps_rate_constant_idx[i_vector];
            rate_constants[i_rate_constant] = rate;
            
            cout << "Rate constant " << rate_constant_names[i_rate_constant] << " set to " << rate << endl;

            // Calculate index for next vector which has to be swept over
            i_parameter = i_parameter/sweeps_rate_constants[i_vector].size();
        }

        // Species numbers
        for(unsigned int i_vector = 0; i_vector < sweeps_species_number.size(); i_vector++){
            
            // Find number of particles and distribute them
            n_particles = sweeps_species_number[i_vector][i_parameter%sweeps_species_number[i_vector].size()];
            i_species = sweeps_species_idx[i_vector];
            initial_species_number[i_species] = n_particles;

            cout << "Species " << species[i_species] << " set to " << n_particles << endl;

            // Calculate index for next vector which has to be swept over
            i_parameter = i_parameter/sweeps_species_number[i_vector].size();
        }

        // Run the model for one parameter set
        run();
    }
}

void Model::define_diffusion_reactions(){

    Reaction* reaction_ptr;
    Event* event_ptr;
    Cell *cell_ptr, *neighbour_ptr;

    int neighbour_idx;

    Cell *phase_sep_neighbour_ptr;
    int phase_sep_neighbour_idx;

    for(int cell_idx = 0; cell_idx < n_cells; cell_idx++){
        cell_ptr = &cells[cell_idx];

        for(unsigned int i_rel_neighbour = 0; i_rel_neighbour < cell_ptr->n_neighbours; i_rel_neighbour++){
            neighbour_idx = cell_ptr->neighbours[i_rel_neighbour];
            neighbour_ptr = &cells[neighbour_idx];

            // Only allow diffusion between cells that have the same dimension. No diffusion between membrane and cytosol
            if(neighbour_ptr->dim != cell_ptr->dim){
                continue;
            }

            for(int i_species = 0; i_species < n_species; i_species++){
                
                // Only allow diffusion of a species that can occur in the cell. Not diffusion of a cytosolic species on a membrane.
                if(species_dimensions[i_species] != cell_ptr->dim){
                    continue;
                }

                event_ptr = &events[n_events];
                reaction_ptr = &diffusion_reactions[i_species];

                event_ptr->reaction_ptr = reaction_ptr;
                event_ptr->target_cells[0] = cell_idx;
                event_ptr->target_cells[1] = neighbour_idx;

                // Check if distance is valid
                if (cell_ptr->neighbour_dist[neighbour_idx] <= 0){
                    cout << "Distance is nonpositive. Distance between cells of same dimensionality must be positive" << endl;
                    cin.get();
                    exit(1);
                }

                // Todo: check if the rate actually exists. The rate constant name might exist, but it may not have a value.
                // I don't know what the basic value is, but it probably isn't pretty.
                event_ptr->rate = rate_constants[reaction_ptr->k_idx]*cell_ptr->neighbour_interface[neighbour_idx]
                                  /(cell_ptr->neighbour_dist[neighbour_idx] * cell_ptr->size);
                event_ptr->base_rate = event_ptr->rate;
                event_ptr->idx = n_events;

                // For selective updating of events. If the cell gets updated, this events also needs to be updated
                // We also need neighbour update, since if cell is at capacity, diffusion from the neighbour cannot occur
                cell_ptr->add_affected_event(event_ptr);
                neighbour_ptr->add_affected_event(event_ptr);

                // If using phase seperation, the (change in) gradient energy is dependent on the concentration of neighbouring cells
                if(phase_separate){

                    // cout << "Adding phase seperation affected events" << endl;
                    
                    // Neighbours of cell
                    for(unsigned int i_phase_sep_neighbour = 0; i_phase_sep_neighbour < cell_ptr->n_neighbours; i_phase_sep_neighbour++){

                        phase_sep_neighbour_idx = cell_ptr->neighbours[i_phase_sep_neighbour];
                        phase_sep_neighbour_ptr = &cells[phase_sep_neighbour_idx];
                        
                        // Don't add to bulk
                        if(phase_sep_neighbour_ptr->dim != cell_ptr->dim){
                            continue;
                        }
                        // Don't add to neighbour twice
                        if(phase_sep_neighbour_ptr != neighbour_ptr){
                            phase_sep_neighbour_ptr->add_affected_event(event_ptr);      
                        }           
                    }

                    // Neighbours of neighbour
                    for(unsigned int i_phase_sep_neighbour = 0; i_phase_sep_neighbour < neighbour_ptr->n_neighbours; i_phase_sep_neighbour++){
                        phase_sep_neighbour_idx = neighbour_ptr->neighbours[i_phase_sep_neighbour];
                        phase_sep_neighbour_ptr = &cells[phase_sep_neighbour_idx];

                        // Don't add to bulk
                        if(phase_sep_neighbour_ptr->dim != neighbour_ptr->dim){
                            continue;
                        }
                        // Don't add to neighbour twice
                        if(phase_sep_neighbour_ptr != cell_ptr){
                            phase_sep_neighbour_ptr->add_affected_event(event_ptr); 
                        }
                    }
                }

                n_events++;
            }
        }
    }


    // cout << "Cell " << cells[0].idx << endl;
    // for(int i_affected_event = 0; i_affected_event < cells[0].n_affected_events; i_affected_event++){
    //     cout << cells[0].affected_events[i_affected_event]->reaction_ptr->reaction_name;
    //     if (cells[0].affected_events[i_affected_event]->reaction_ptr->is_diffusion_like){
    //         cout << "From cell " << cells[0].affected_events[i_affected_event]->target_cells[0] << " to " << cells[0].affected_events[i_affected_event]->target_cells[1] << endl;
    //     } else {cout << endl;}

    //     if(i_affected_event % 100 == 0){
    //         cin.get();
    //     }
    // }
    // cout << "end" << endl;
    // cin.get();

}

void Model::define_normal_reactions(){

    Cell* cell_ptr;
    Reaction* reaction_ptr;
    Event* event_ptr;
    
    for(int i_reaction = 0; i_reaction < n_reactions; i_reaction++){

        // Get data from reaction class for species and rate
        reaction_ptr = &reactions[i_reaction];   

        for(int cell_idx = 0; cell_idx < n_cells; cell_idx++){

            cell_ptr = &cells[cell_idx];

            event_ptr = &events[n_events];
            event_ptr->reaction_ptr = reaction_ptr;

            for(int i_reaction_species = 0; i_reaction_species < reaction_ptr->n_reaction_species; i_reaction_species++){
                event_ptr->target_cells[i_reaction_species] = cell_idx;
            }

            // Set base event rate
            event_ptr->rate = rate_constants[reaction_ptr->k_idx];

            // Correct for cell size
            event_ptr->rate *= cells[cell_idx].size * free_size_fraction[cells[cell_idx].dim];
            // cout << "Free size fraction: " << free_size_fraction[cells[cell_idx].dim] << endl;

            for(int i_reactant = 0; i_reactant < reaction_ptr->n_reactants; i_reactant++){
                event_ptr->rate /= cells[cell_idx].size * free_size_fraction[cells[cell_idx].dim];
                // cout << "Free size fraction: " << free_size_fraction[cells[cell_idx].dim] << endl;
            }

            event_ptr->base_rate = event_ptr->rate;
            event_ptr->idx = n_events;
            cell_ptr->add_affected_event(event_ptr);

            // Update event count
            n_events++;
        }
    }

}

void Model::define_boundary_reactions(){

    Reaction* reaction_ptr;
    Event* event_ptr;
    Cell* cell_ptr, *neighbour_ptr;

    int neighbour_idx;
    int species_idx, species_dim;

    for(int i_reaction = 0; i_reaction < n_boundary_reactions; i_reaction++){

        // Get data from reaction class for species and rate
        
        reaction_ptr = &boundary_reactions[i_reaction];   
        

        // We assume that each boundary reaction only involves 2 cells. If we go over all unidirectional edges,
        // and check if we can assign reactants to the 2 cells, we add each reaction twice
        // E.g. membrane dissociation. First we find a 2D cell with 3D neighbour, and add reaction appropriately.
        // Then we find 3D cell with samen 2D neighbour, and again add reaction, but this is the same reaction!
        // Hence, we only check if the dimension of the first reactant fits.

        // Go over all neighbours. 
        for(int cell_idx = 0; cell_idx < n_cells; cell_idx++){
            cell_ptr = &cells[cell_idx];

            // First reactant dimension does not match. Continue so there is no double counting.
            if (species_dimensions[reaction_ptr->reaction_species[0]] != cell_ptr->dim){
                continue;
            }

            for(unsigned int i_neighbour = 0; i_neighbour < cell_ptr->n_neighbours; i_neighbour++){

                neighbour_idx = cell_ptr->neighbours[i_neighbour];
                neighbour_ptr = &cells[neighbour_idx];

                // Neighbour must be a different type
                if (cell_ptr->dim == neighbour_ptr->dim){
                    continue;
                }

                event_ptr = &events[n_events];
                event_ptr->reaction_ptr = reaction_ptr;

                // Setting target cells
                for (int i_reaction_species = 0; i_reaction_species < reaction_ptr->n_reaction_species; i_reaction_species++){
                    if (species_dimensions[reaction_ptr->reaction_species[i_reaction_species]] == cell_ptr->dim){
                        event_ptr->target_cells[i_reaction_species] = cell_idx;
                    } else if (species_dimensions[reaction_ptr->reaction_species[i_reaction_species]] == neighbour_ptr->dim){
                        event_ptr->target_cells[i_reaction_species] = neighbour_idx;  
                    } else{
                        cout << "The dimension of reaction species " << species[reaction_ptr->reaction_species[i_reaction_species]] << " does not equal the dimension of either the original or neighbour cell" << endl;
                        cout << "This can happen if there are more than 2 different dimensions, but if there are only 2 different dimensions, something went wrong" << endl;
                        cin.get();
                    }
                }

                // TODO add in more documentation for rate correction.
                // When adding reaction, need to add which right hand reagent k is relative to.
                // In file, denote this by a new line.
                
                // Set base event rate
                event_ptr->rate = rate_constants[reaction_ptr->k_idx];

                // Correct for volume or area
                // All constants are relative to the membrane species.
                if (cell_ptr->dim < neighbour_ptr->dim){
                    event_ptr->rate *= cell_ptr->size * free_size_fraction[cell_ptr->dim];
                } else if (neighbour_ptr->dim < cell_ptr->dim){
                    event_ptr->rate *= neighbour_ptr->size * free_size_fraction[neighbour_ptr->dim];
                } else{
                    cout << "Cell " << cell_ptr->idx << ". Neighbour " << neighbour_ptr->idx << endl;
                    cout << "Dimesions " << cell_ptr->dim << " and " << neighbour_ptr->dim << endl;
                    cout << "Cell and neighbour dimensions are the same, should not be possible" << endl;
                    cin.get();
                    exit(1);
                }

                bool includes_bulk_component = false;

                for(int i_reactant = 0; i_reactant < reaction_ptr->n_reactants; i_reactant++){
                    species_idx = reaction_ptr->reaction_species[reaction_ptr->reactants[i_reactant]];
                    species_dim = species_dimensions[species_idx];

                    // Temporary code
                    // If event includes a bulk component, it needs to be added to affected events of the bulk
                    // Else, only the the affected events of the other cell.
                    if (species_dim == 3){
                        includes_bulk_component = true;
                    }
                    
                    // Dividing by cell size for each reactant
                    if(species_dim == cell_ptr->dim) {
                        event_ptr->rate /= cell_ptr->size * free_size_fraction[cell_ptr->dim];
                        // cout << "Free size fraction: " << free_size_fraction[cell_ptr->dim] << endl;
                    } else if(species_dim == neighbour_ptr->dim) {
                        event_ptr->rate /= neighbour_ptr->size * free_size_fraction[neighbour_ptr->dim];
                        // cout << "Free size fraction: " << free_size_fraction[neighbour_ptr->dim] << endl;
                    } else {
                        cout << "The species dimension " << species_dim << "is not equal to either of the cell dimensions. Check the input file." << endl;
                        cin.get();
                        exit(1);
                    }
                }

                event_ptr->base_rate = event_ptr->rate;
                event_ptr->idx = n_events;
                
                // Temporary code
                // If reaction includes bulk component, add to affected event bulk for propensity
                // and to other cell for checking capacity
                if (includes_bulk_component){
                    cell_ptr->add_affected_event(event_ptr);
                    neighbour_ptr->add_affected_event(event_ptr);
                } else{
                    if (cell_ptr->dim == 3 && neighbour_ptr->dim == 2){
                        neighbour_ptr->add_affected_event(event_ptr);
                    } else if (neighbour_ptr->dim == 3  && cell_ptr->dim == 2){
                        cell_ptr->add_affected_event(event_ptr);
                    } else{
                        cout << "Neither cell nor neighbour is bulk" << endl;
                        cin.get(); exit(1);
                    }
                }

                // Update event count
                n_events++;                
            }
        }
    }
}


Vertex::Vertex(){}
Vertex::Vertex(double x, double y, double z){
    this->x = x;
    this->y = y;
    this->z = z;
}

double Vertex::distance(const Vertex &neighbour){

    return sqrt(pow(x - neighbour.x, 2) + pow(y - neighbour.y, 2) + pow(z - neighbour.z, 2));
}

double Vertex::abs(){

    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

Vertex Vertex::operator+(const Vertex &v2){
    Vertex v3(this->x + v2.x, this->y + v2.y, this->z + v2.z);
    return v3;
}

Vertex Vertex::operator-(const Vertex &v2){
    Vertex v3(this->x - v2.x, this->y - v2.y, this->z - v2.z);
    return v3;
}

Vertex& Vertex::operator/=(const float &div){
    this->x /= div;
    this->y /= div;
    this->z /= div;
    return *this;
}

Vertex Vertex::operator/(const float &div){
    Vertex v3 = *this;
    return v3 /= div;
}

Vertex Vertex::cross(const Vertex &v2){
    Vertex v3(this->y*v2.z - this->z*v2.y,
              this->z*v2.x - this->x*v2.z,
              this->x*v2.y - this->y*v2.x);
    return v3;
}

double Vertex::dot(const Vertex &v2){
    return this->x*v2.x + this->y*v2.y + this->z*v2.z;
}

Edge::Edge(Cell &c1, Cell &c2){
    this->c1 = &c1;
    this->c2 = &c2;
}

Face::Face(Cell *c1, Cell *c2, Cell *c3){
    this->cells.push_back(c1);
    this->cells.push_back(c2);
    this->cells.push_back(c3);
}

void Face::generate_angles_and_lengths(){

    // cout << "Triangle " << cells[0]->idx << ", " << cells[1]->idx << ", " << cells[2]->idx << endl;
    // cout << "Positions cell 1: " << cells[0]->pos.x << ", " << cells[0]->pos.y << ", " << cells[0]->pos.z << endl;
    // cout << "Positions cell 2: " << cells[1]->pos.x << ", " << cells[1]->pos.y << ", " << cells[1]->pos.z << endl;
    // cout << "Positions cell 3: " << cells[2]->pos.x << ", " << cells[2]->pos.y << ", " << cells[2]->pos.z << endl;

    lengths.reserve(cells.size());
    theta_terms.reserve(cells.size());


    for(unsigned int i_cell = 0; i_cell < cells.size(); i_cell++){
        lengths.push_back((cells[(i_cell+1)%3]->pos).distance(cells[(i_cell+2)%3]->pos));
        theta_terms.push_back((cells[(i_cell+1)%3]->pos - cells[i_cell]->pos).dot(cells[(i_cell+2)%3]->pos - cells[i_cell]->pos));
    }

    area = (cells[1]->pos - cells[0]->pos).cross(cells[2]->pos - cells[0]->pos).abs()/2;

    // for(unsigned int i_cell = 0; i_cell < cells.size(); i_cell++){
    //     cout << "Length: " << lengths[i_cell] << endl;
    //     cout << "Area: " << area << endl;
    //     cout << "Expected area: " << (cells[(i_cell+1)%3]->pos - cells[i_cell]->pos).cross(cells[(i_cell+2)%3]->pos - cells[i_cell]->pos).abs()/2 << endl;
    //     cout << "Theta term: " << theta_terms[i_cell] << endl;
    //     cout << "Expected theta term: " << (pow(lengths[(i_cell+1)%3], 2) + pow(lengths[(i_cell+2)%3], 2) - pow(lengths[i_cell], 2))/2 << endl;
    // }
    // cin.get();

}

// void Face::generate_edge_normals(){

//     // Finding normal of entire face
//     Vertex face_normal = (cells[1]->pos - cells[0]->pos).cross(cells[2]->pos - cells[1]->pos);

//     edge_normals.reserve(cells.size());
    
//     // Finding normal of edge, always pointing outward from face.
//     for(unsigned int i_cell = 0; i_cell < cells.size(); i_cell++){
//         Vertex edge_normal = (cells[i_cell]->pos - cells[(i_cell + 1)%cells.size()]->pos).cross(face_normal);
//         edge_normal /= edge_normal.abs();
//         edge_normals.push_back(edge_normal);
//         // cout << "Edge normal: " << edge_normal.x << ", " << edge_normal.y << "," << edge_normal.z << endl;
//         // cout << "Size: " << edge_normal.abs() << endl;
//         // cout << "Dot: " << (edge_normal.dot(cells[i_cell]->pos - cells[(i_cell + 1)%cells.size()]->pos)) << endl;

//     }
// }