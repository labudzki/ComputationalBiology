#include "SpatialSSA_model.hpp"

// Only if not on cluster
// #include "SpatialSSA_rendering.hpp"

void run(Model* model_ptr){

    // In case parameter sweep is added through input file
    // model_ptr->do_parameter_sweep();

    // In case parameter sweep is done through cluster
    model_ptr->run();

    cout << "Press any key to continue";
    cin.get();
}

int main(int argc, char** argv) {

    // Loading in directory name
    string dir_name = "";
    if(argc > 1){
        dir_name = argv[1];
    } else{
        dir_name = "Data";
    }


    // Creating the model. Put it in heap memory, so that it can be larger. 
    Model* model_ptr = new Model;
    model_ptr->dir_name = dir_name;
    cout << "Directory name: " << model_ptr->dir_name << endl;

    // Loading in sweeps
    if(argc > 2 && argc%2 == 1){
        cout << "There are unmatched variables in the input" << endl;
        cin.get(); exit(1);
    }

    // Counting number of variables
    int n_variables;
    if (argc > 2){n_variables = (argc - 2)/2;}
    else {n_variables = 0;}

    vector<string> variable_names; variable_names.reserve(n_variables);
    vector<double> variable_vals; variable_vals.reserve(n_variables);
    
    // Storing sweep variables
    for(int i_variable = 0; i_variable < n_variables; i_variable++){
        cout << "Saving variable: " << argv[2*i_variable + 2] << endl;
        cout << "With value: " << stod(argv[2*i_variable + 3]) << endl;
        variable_names.push_back(argv[2*i_variable + 2]);
        variable_vals.push_back(stod(argv[2*i_variable + 3]));
        cout << "Variable saved" << endl;
    }

    

    // ----------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------Phase seperation setup---------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------

    // // These work well for phase seperation
    // // model_ptr->create_spherical_grid(162, 0.282);
    // // model_ptr->create_spherical_grid(642, 0.282);
    // model_ptr->create_spherical_grid(2562, 0.282);
    // // model_ptr->create_cartesian_grid(1, 1, 1, 50, 50, true);

    // // model_ptr->load_model("./Models/diffusion_absorption_model.txt");
    // // model_ptr->load_initial_conditions("./Models/diffusion_absorption_setup.txt");
    // model_ptr->load_model("./Models/diffusion_model.txt");
    // model_ptr->load_initial_conditions("./Models/diffusion_setup.txt");


    // model_ptr->phase_separate = true;

    // model_ptr->verbose = false;
    // model_ptr->timing = true;

    // model_ptr->record_indv_cells = false;
    // model_ptr->record_propensities = false;

    // model_ptr->break_iter = 5000;
    // model_ptr->n_iterations = 1000000;

    // // Without indv cell and propensities 300ms for 10 records on 642 cells
    // // With individual cells 7200ms for 10 records on 642 cells
    // // With propensities 1050ms for 10 records on 642 cells

    // // 100 million iterations, break every 100 thousand.
    // // model_ptr->n_iterations = 10000000;


    // ----------------------------------------------------------------------------------------------------------------
    // --------------------------------------------Yeast polarization setup--------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------

    // Replicates work by Frey et al.
    // model_ptr->create_spherical_grid(642, 3.5);
    model_ptr->create_spherical_grid(2562, 3.5);

    model_ptr->load_model("./Models/yeast_model.txt");
    model_ptr->load_initial_conditions("./Models/yeast_setup.txt");

    // Model 1 inactivated. No GAP saturation by high GAP activity
    // Should polarize. Does polarize.
    // model_ptr->load_model("./Models/yeast_model_1_GAP_activity.txt");
    // model_ptr->load_initial_conditions("./Models/yeast_setup_1_GAP_activity.txt");

    // Model 2 inactivated. No Cdc42 recruitment by Bem1 by Bem1 knockout
    // Should not polarize in this regime. Polarizes a little, maybe.
    // model_ptr->load_model("./Models/yeast_model_2_Bem1-knockout.txt");
    // model_ptr->load_initial_conditions("./Models/yeast_setup_2_Bem1-knockout.txt");

    // Model 3 inactivated. No Cdc42 movement, since it is immobile and already distributed on the membrane.
    // Should polarize. Does not polarize, not noticeable. Maybe very little.
    // Does not polarize if starting all on membrane. Does polarize when starting in bulk
    // Current model starts all on membrane!
    // model_ptr->load_model("./Models/yeast_model_3_Cdc42-immobile.txt");
    // model_ptr->load_initial_conditions("./Models/yeast_setup_3_Cdc42-immobile.txt");


    model_ptr->phase_separate = false;

    model_ptr->verbose = false;
    model_ptr->timing = true;

    model_ptr->record_indv_cells = false;
    model_ptr->record_propensities = false;

    // Simulate 10 seconds, record every 5000 iterations.
    model_ptr->t_max        = 200;
    model_ptr->break_iter   = 5000;
    model_ptr->n_iterations = 2000000; // 2 million for 2562
    // model_ptr->n_iterations = 500000; // 500 thousand for 642

    model_ptr->n_repeats = 1;

    // ===================================
    // Particles per cell is still strange
    // ===================================
    

    // // TODO: put timing variables back in correct place.
    // time(&t_start);
    // time(&t_end);


    // Setting variables from input line
    int species_idx;

    for(int i_variable = 0; i_variable < n_variables; i_variable++){

        cout << "Variable name: " << variable_names[i_variable] << endl;
        cout << "Variable value: " << variable_vals[i_variable] << endl;

        // Go over rate constants
        if(model_ptr->rate_constant_map.count(variable_names[i_variable])){
            // If name occurs in rate constant map names.
            model_ptr->rate_constants[model_ptr->rate_constant_map[variable_names[i_variable]]] = variable_vals[i_variable];
            break;
        } 
        else {
            // If not a already define rate constant, it must be a species. If not, throw error (get_species_idx does this)
            species_idx = model_ptr->get_species_idx(variable_names[i_variable]);
            model_ptr->initial_species_number[species_idx] = int(variable_vals[i_variable]);
        }
    }

    
    // ===================================
    //           Rendering mode
    // ===================================

    // Also change imports
    // -------------------

    // sf::Thread thread(&run, model_ptr); //Mac does not allow the opening of windows that are not in the main thread, for ... reasons, so the simulations occurs in a separate thread
    // thread.launch();
    // // render(model_ptr, false, true);
    // render(model_ptr, true, false);
    
    // ===================================
    //         Computing mode
    // ===================================

    run(model_ptr);
    
    return 0;
}

// model_ptr->cells[0].set_species_number(0, 10000);
// model_ptr->cells[0].set_species_number(1, 130);
// model_ptr->cells[1].set_species_number(0, 500);

// TODO: also check events in event creation loop.
// Check events with correct rate constants, targets, everything!!!
// Checking reactions

// Reaction* reaction_ptr;
// for(int i = 0; i < model_ptr->n_reactions; i++){
//     reaction_ptr = &(model_ptr->reactions[i]);

//     cout << reaction_ptr->reaction_name << endl;

//     for(int j = 0; j < reaction_ptr->n_reaction_species;j++){
//         cout << "Species "<< reaction_ptr->species[j] << " changes by " <<
//         reaction_ptr->state_change[j] << " with rate " << reaction_ptr->k << endl;
//     }
// }

// Checking Reactions
// for(int j = 0; j<3; j++){
// for(int i = 0; i<model.reactions[j].n_reaction_species; i++){
//     cout << "species "<< model.reactions[j].species[i] << " changes by " << model.reactions[j].state_change[i] << "\n";
// }
// }

// Checking neighbours
// for(int i = 0; i<model_ptr->n_cells; i++){

//     cout << "cell " << i << " neighbours \n";

//     for(int j = 0; j<model_ptr->cells[i].n_neighbours; j++){
        
//         cout << "Neighbour: " << model_ptr->cells[i].neighbours[j] << endl;
//         cout << "Distance: " << model_ptr->cells[i].neighbour_dist[j] << endl;

//     }
//     cout << endl;
// }

