#include "SpatialSSA_model.hpp"

void check_valid_name(string name){
    
    const string forbidden_chars = "<>:\"/\\|?*,.";
    for(unsigned int i_char = 0; i_char < forbidden_chars.length(); i_char++){
        if (name.find(forbidden_chars[i_char]) != string::npos) {
            cout << "Name " << name << " contains forbidden character: " << forbidden_chars[i_char] << endl;
            cin.get();
            exit(1);
        }
    }
}

void Model::load_model(const string file){
    
    ifstream setup_file;
    setup_file.open(file);
    if (!setup_file){
        cout << "File " << file << " does not exist" << endl;
        cin.get();
        exit(1);
    }

    // Variables for data parsing
    stringstream stream;
    string line;

    string species_name, rate_constant_name, reaction_name;
    
    int species_dimensions;
    int indv_change;
    int n_reaction_species;
   
    string reactants[MAX_REACTION_SIZE];
    int state_change[MAX_REACTION_SIZE];


    // Loading in model name
    while (!setup_file.eof()) {

        getline(setup_file, line);

        // If comment or empty line, skip line
        if(line[0] == '/' || line.length() == 0){
            continue;
        }

        // If symbol indicates next section, break this loop
        if(line[0] == '>'){
            break;
        }

        model_name = line;
    }

    // Loading in rate constants
    while (!setup_file.eof()) {

        // Save the line in "line".
        getline(setup_file, line); 

        // If comment or empty line, skip line
        if(line[0] == '/' || line.length() == 0){
            continue;
        }

        // If symbol indicates next section, break this loop
        if(line[0] == '>'){
            break;
        }

        // Tell stream that it has not reached the end.
        stream.clear();

        stream << line;
        stream >> rate_constant_name;
        cout << "Reading rate constant " << rate_constant_name;
        cout << endl;

        // Sanitizing species name:
        check_valid_name(rate_constant_name);
        add_rate_constant(rate_constant_name);
    }

    cout << "Number of rate constants: " << n_rate_constants << endl;

    // Loading in species data
    while (!setup_file.eof()) {

        // Save the line in "line".
        getline(setup_file, line); 

        // If comment or empty line, skip line
        if(line[0] == '/' || line.length() == 0){
            continue;
        }

        // If symbol indicates next section, break this loop
        if(line[0] == '>'){
            break;
        }

        // Tell stream that it has not reached the end.
        stream.clear();

        stream << line;
        stream >> species_name >> species_dimensions >> rate_constant_name;
        cout << "Reading species " << species_name << ", dimension: " << species_dimensions << ", diffusion constant: " << rate_constant_name << endl;

        // Sanitizing species name:
        check_valid_name(species_name);
        add_species(species_name, species_dimensions, rate_constant_name);
    }

    // Loading in reaction data
    while (!setup_file.eof()) {

        // Save the line in "line".
        getline(setup_file, line); 

        // If comment or empty line, skip line
        if(line[0] == '/' || line.length() == 0){
            continue;
        }

        // Get reaction name
        reaction_name = line;
        cout << "Reading reaction " << reaction_name << endl;

        getline(setup_file, line); 

        stream.clear();

        // Get reactants
        stream << line;

        n_reaction_species = 0;
        while (!stream.eof()){

            if(n_reaction_species == MAX_REACTION_SIZE){
                cout << "Number of reactants in reaction " << reaction_name << " exceeds maximum reaction size" << endl;
                cin.get();
                exit(1);
            }

            stream >> species_name >> indv_change;

            cout << "Species " << species_name << " , change: " << indv_change << endl;
            reactants[n_reaction_species] = species_name;
            state_change[n_reaction_species] = indv_change;
            n_reaction_species++;
        }

        // Get reaction rate
        stream.clear();
        getline(setup_file, line);
        stream << line;
        stream >> rate_constant_name;
        cout << "Rate constant: " << rate_constant_name << endl;
        add_reaction(reactants, state_change, rate_constant_name, n_reaction_species, reaction_name);


        // TODO: add check that some combination of reactants remain constant.
    }
}

void Model::load_initial_conditions(const string file){

    // TODO IDEA. Allow marking of cells with flags. Then add species randomly or uniformly over those flags.
    // TODO: use DIR from old code

    ifstream setup_file;
    setup_file.open(file);
    if (!setup_file){
        cout << "File " << file << " does not exist" << endl;
        cin.get();
        exit(1);
    }

    stringstream stream;
    string line;
    string species_name, rate_constant_name;
    string sweep_type;

    // For initial particle distribution loading
    // For rate and diffusion constant loading
    int i_species, i_rate_constant;
    int n_start, n_end, n_points;
    double k_start, k_end;


    cout << "Project name: " << endl;

    // Loading the project name to be used in save files.
    while (!setup_file.eof()) {

        getline(setup_file, line);

        // If comment or empty line, skip line
        if(line[0] == '/' || line.length() == 0){
            continue;
        }

        // If symbol indicates next section, break this loop
        if(line[0] == '>'){
            break;
        }

        run_name = line;
        cout << run_name << endl;
        cout << endl;

    }

    check_valid_name(run_name);

    cout << "Initial conditions: " << endl;

    // Loading the initial number of particles.
    while (!setup_file.eof()) {

        getline(setup_file, line);

        // If comment or empty line, skip line
        if(line[0] == '/' || line.length() == 0){
            continue;
        }

        // If symbol indicates next section, break this loop
        if(line[0] == '>'){
            break;
        }

        // Resetting stream and variables
        species_name = sweep_type =  "";
        n_start = n_end = n_points = 0;
        stream.clear();

        // Adding line to stream and extracting data.
        stream << line;
        stream >> species_name >> n_start >> n_end >> sweep_type >> n_points;

        // Incorrect syntax
        if(n_end != 0 && sweep_type == ""){
            cout << "For " << species_name << " add desired type of parameter sweep (lin or log) to line" << endl;
            cin.get();
            exit(1);
        }

        // Negative particle numbers
        if(n_start < 0 || n_end < 0){
            cout << "For " << species_name << ", particle number may not be negative" << endl;
            cin.get();
            exit(1);
        }

        // No sweep, just distribute particles.
        if(sweep_type == ""){
            i_species = get_species_idx(species_name);
            initial_species_number[i_species] = n_start;
            cout << species_name << ": " << n_start << endl;
        }

        // Particle number sweep
        else{

            // Safety check on particle number
            if(n_start >= n_end){
                cout << "For " << species_name << ", starting number of particles (" << n_start << ") must be smaller than final number of particles (" << n_end << ")" << endl;
                cin.get();
                exit(1);
            }
            if(sweep_type == "lin" || sweep_type == "log"){
                i_species = get_species_idx(species_name);
                cout << species_name << ": " << sweep_type << " sweep from " << n_start << " to " << n_end << endl;
                add_species_sweep(i_species, n_start, n_end, n_points, sweep_type);

            } else {
                cout << "For " << species_name << ", sweep type " << sweep_type << " is not recognised" << endl;
                cin.get();
                exit(1);
            }
        }
    }
    cout << endl;

    // Load rate constants values
    cout << "Rate constants and diffusion coefficients:" << endl;

    while (!setup_file.eof()) {

        getline(setup_file, line);

        // If comment or empty line, skip line
        if(line[0] == '/' || line.length() == 0){
            continue;
        }

        // If symbol indicates next section, break this loop
        if(line[0] == '>'){
            break;
        }

        // Resetting stream and variables
        rate_constant_name = sweep_type =  "";
        k_start = k_end = 0;
        n_points = 0;
        stream.clear();

        // Pushing line to stream and extracting name and value
        stream << line;
        stream >> rate_constant_name >> k_start >> k_end >> sweep_type >> n_points;


        // Incorrect syntax
        if(k_end != 0 && sweep_type == ""){
            cout << "For " << rate_constant_name << " add desired type of parameter sweep (lin or log) to line" << endl;
            cin.get();
            exit(1);
        }

        // Negative rate constant
        if(k_start < 0 || k_end < 0){
            cout << "For " << rate_constant_name << ", value may not be negative" << endl;
            cin.get();
            exit(1);
        }

        // No sweep, just set rate constant.
        if(sweep_type == ""){
            i_rate_constant = get_rate_idx(rate_constant_name);
            cout << rate_constant_name << " = " << k_start << endl;
            rate_constants[i_rate_constant] = k_start;
        }

        // Rate constant sweep
        else{

            // Safety check on particle number
            if(k_start > k_end){
                cout << "For " << rate_constant_name << ", starting value (" << k_start << ") must be smaller than final value (" << k_end << ")" << endl;
                cin.get();
                exit(1);
            }
            if(sweep_type == "lin" || sweep_type == "log"){
                i_rate_constant = get_rate_idx(rate_constant_name);
                cout << rate_constant_name << ": " << sweep_type << " sweep from " << k_start << " to " << k_end << endl;
                add_rate_constant_sweep(i_rate_constant, k_start, k_end, n_points, sweep_type);

            } else {
                cout << "For " << rate_constant_name << ", sweep type " << sweep_type << " is not recognised" << endl;
                cin.get();
                exit(1);
            }
        }
    }

    cout << endl;
}

// Recording

void Model::record_concentrations(const int n_iter) {

    // Recording all species in a cell, for all cells
    if(record_indv_cells){
        for(int i_cell = 0; i_cell < n_cells; i_cell++){
            record_all_species_in_cell(n_iter, i_cell);
        }
    }

    // Recording species distribution over all cells, for all species.
    if(record_indv_species){
        for(int i_species = 0; i_species < n_species; i_species++){
            record_all_cells_for_species(i_species, n_iter);
        }
    }

    if(record_reactions){
        record_all_reactions(n_iter);
    }

    if(record_propensities){
        record_all_propensities(n_iter);
    }

    record_mu(n_iter);

    initial_recording = false;

}


void Model::record_all_species_in_cell(const int n_iter, const int i_cell) {
	
    // Recording the number of components as a function of time in a given cell.
    
    // Setting file name to include project and cell index.
    string file_name, cell_name;
    cell_name = "_components_in_cell_" + to_string(i_cell);
    file_name = project_name + cell_name + ".csv";

	ofstream datafile;

	if(initial_recording) {
		
		datafile.open(file_name);
		// datafile << "#Evolution of concentrations of " << model_name << " species in one grid cell" << endl;
		// datafile << "#The components on grid cell " << i_cell << endl << endl;

		datafile << setw(20) << "Step (n)," << setw(19) << "Time (t)";

        for(int i_species = 0; i_species < n_species; i_species++){
            datafile  << "," << setw(19) << species[i_species];
        }

		datafile << "\n";
				
	}
	
	else {datafile.open(file_name, ios::app);}
	
	datafile << setw(19) << n_iter << "," << setw(19) << t;
	for (int i_species = 0; i_species < n_species; i_species++) {
		datafile  << "," << setw(19) << cells[i_cell].state[i_species];
	}

	datafile << "\n";
	datafile.close();
}

void Model::record_all_cells_for_species(const int i_species, const int n_iter) {
	
    // Recording the number of components as a function of time in a given cell.
    
    // Setting file name to include project and cell index.
    string file_name, species_name, species_description;
    int species_distribution[MAX_CELLS];

    species_name = species[i_species];
    species_description = "_concentration_of_species_" + species_name;
    file_name = project_name + species_description + ".csv";


    for(int i_cell = 0; i_cell < n_cells; i_cell++){
        species_distribution[i_cell] = cells[i_cell].state[i_species];
    }

	ofstream datafile;

	if(initial_recording) {
		
		datafile.open(file_name);
		// datafile << "#Evolution of concentrations of " << species_name << " in " << model_name << " in all cells" << endl;

		datafile << "Step (n),Time (t)";

        for(int i_cell = 0; i_cell < n_cells; i_cell++){
            datafile  << ",Cell " + to_string(i_cell);
        }
		datafile << "\n";
				
	}
	
	else {datafile.open(file_name, ios::app);}
	
	datafile << n_iter << "," << t;
	for (int i_cell = 0; i_cell < n_cells; i_cell++) {
		datafile  << "," << species_distribution[i_cell];
	}

	datafile << '\n';
	datafile.close();
}

void Model::record_all_propensities(const int n_iter){

    string file_name;
    file_name = project_name + "_propensities.csv";

    ofstream datafile;

    if(initial_recording) {
		
		datafile.open(file_name);
		// datafile << "#Number of reactions of each type in " << model_name << endl;
		datafile << setw (20) << "Step (n)," << setw (19) << "Time (t)";

        for(int i_event = 0; i_event < n_events; i_event++){
            datafile << "," << setw(19) << "Event " + to_string(i_event);
        }
		datafile << "\n";	
	}

    else {datafile.open(file_name, ios::app);}

	datafile << setw (19) << n_iter << "," << setw (19) << t;
    for(int i_event = 0; i_event < n_events; i_event++){
        datafile << "," << setw(19) << propensities[i_event];
    }

    datafile << endl;
	datafile.close();
}


void Model::record_all_reactions(const int n_iter){

    string file_name;
    file_name = project_name + "_reaction_counts.csv";

    ofstream datafile;

    if(initial_recording) {
		
		datafile.open(file_name);
		// datafile << "#Number of reactions of each type in " << model_name << endl;
		datafile << setw (20) << "Step (n)," << setw (19) << "Time (t)";

        for(int i_reaction = 0; i_reaction < n_reactions; i_reaction++){
            datafile << "," << setw(49) << reactions[i_reaction].reaction_name;
        }

        for(int i_boundary_reaction = 0; i_boundary_reaction < n_boundary_reactions; i_boundary_reaction++){
            datafile << "," << setw(49) << boundary_reactions[i_boundary_reaction].reaction_name;
        }

        for(int i_species = 0; i_species < n_species; i_species++){
            datafile << "," << setw(49) << diffusion_reactions[i_species].reaction_name;
        }
		datafile << "\n";	
	}

    else {datafile.open(file_name, ios::app);}

	datafile << setw (19) << n_iter << "," << setw (19) << t;
    for(int i_reaction = 0; i_reaction < n_reactions; i_reaction++){
        datafile << "," << setw(49) << reactions[i_reaction].count;
    }
    for(int i_boundary_reaction = 0; i_boundary_reaction < n_boundary_reactions; i_boundary_reaction++){
        datafile << "," << setw(49) << boundary_reactions[i_boundary_reaction].count;
    }
    for(int i_species = 0; i_species < n_species; i_species++){
        datafile << "," << setw(49) << diffusion_reactions[i_species].count;
    }

    datafile << endl;
	datafile.close();
}


void Model::record_mu(const int n_iter){

    string file_name;
    file_name = project_name + "_mu.txt";

    ofstream datafile;

    if(initial_recording) {
		
		datafile.open(file_name);
		datafile << setw (19) << "Step (n)";

        for(int i = 0; i < break_iter; i++){
            datafile << "," << setw(8) << to_string(i);
        }
		datafile << "\n";	
	}

    else {
        datafile.open(file_name, ios::app);

        datafile << setw (19) << n_iter;
        for(int i = 0; i < break_iter; i++){
            datafile << "," << setw(8) << mu_history[i];
        }
    }

    mu_history.clear();
    mu_history.reserve(break_iter);

    datafile << endl;
	datafile.close();

}


void Model::record_event_data(){
    // Recording the target cells and state changes of each event

    Event *event_ptr;
    Reaction *reaction_ptr;

    string file_name;
    file_name = project_name + "_event_data.txt";

    ofstream datafile;
    datafile.open(file_name);

    datafile << setw(8) << "Event," << setw(50) << "Reaction name," << setw(20) << "Number of species," << setw(32) << "(Cell; species; state change)," << "\n";
    for(int i_event = 0; i_event < n_events; i_event++){
        event_ptr = &events[i_event];
        reaction_ptr = event_ptr->reaction_ptr;
        datafile << setw(8) << to_string(event_ptr->idx) + "," << setw(50) << reaction_ptr->reaction_name + "," << setw(20) << to_string(reaction_ptr->n_reaction_species) + ',';
        for(int i_reaction_species = 0; i_reaction_species < reaction_ptr->n_reaction_species; i_reaction_species++){
            datafile << setw(32) << to_string(event_ptr->target_cells[i_reaction_species]) + "; " + to_string(reaction_ptr->reaction_species[i_reaction_species]) + "; "
                                     + to_string(reaction_ptr->state_change[i_reaction_species]) + ",";
        }
        datafile << "\n";
    }

    datafile << endl;
	datafile.close();



}

void Model::record_cell_data(){

    Cell *cell_ptr;

    string file_name;
    file_name = project_name + "_cell_properties.txt";

    ofstream datafile;
    datafile.open(file_name);

    // Logging neighbour independent data
    datafile << setw(5) << "Cell," << setw(12) << "Dimension," << setw(14) << "Size," << setw(40) << "Coordinates" << "\n";
    for(int i_cell = 0; i_cell < n_cells; i_cell++){
        datafile << setw(5) << to_string(cells[i_cell].idx) + "," << setw(12) << to_string(cells[i_cell].dim) + ","
                 << setw(14) << to_string(cells[i_cell].size) + ","
                 << setw(40) << to_string(cells[i_cell].pos.x) + "; " + to_string(cells[i_cell].pos.y) + "; " + to_string(cells[i_cell].pos.z) << "\n";
    }
    datafile << "\n";

    // Logging neighbour data
    datafile << setw(5) << "Cell," << setw(25) << "Number of neighbours," << setw(30) << "(Neighbour; dist; int)," << "\n";
    for(int i_cell = 0; i_cell < n_cells; i_cell++){
        cell_ptr = &cells[i_cell];
        datafile << setw(5) << to_string(cell_ptr->idx) + "," << setw(24) << to_string(cell_ptr->n_neighbours);
        for(unsigned int i_neighbour = 0; i_neighbour < cell_ptr->n_neighbours; i_neighbour++){
            datafile << "," << setw(29) << to_string(cell_ptr->neighbours[i_neighbour]) + "; " + to_string(cell_ptr->neighbour_dist[cell_ptr->neighbours[i_neighbour]]) + "; "
                                     + to_string(cell_ptr->neighbour_interface[cell_ptr->neighbours[i_neighbour]]);
        }
        datafile << "\n";
    }

    datafile << endl;
	datafile.close();

}

void Model::record_parameters(){

    string file_name;
    file_name = project_name + "_parameter_values.txt";

    ofstream datafile;
    datafile.open(file_name);

    datafile << "Values of rate constants and diffusion coefficients in " << model_name << " run " << i_run << endl;
    for(int i_parameter = 0; i_parameter < n_rate_constants; i_parameter++){
        datafile << setw(8) << left << rate_constant_names[i_parameter];
        datafile << setw(5) << rate_constants[i_parameter] << endl;
    }

    datafile << endl;

    datafile << "Values for initial number of species" << endl;
    for(int i_species = 0; i_species < n_species; i_species++){
        datafile << setw(15) << left << species[i_species];
        datafile << setw(5) << initial_species_number[i_species] << endl;
    }

	datafile << endl;
	datafile.close();
}


void Model::record_initial_conditions(){

    string file_name;
    file_name = project_name + "_initial_conditions.csv";

    ofstream datafile;
    datafile.open(file_name);

    datafile << "Cell";
    for (int i_species = 0; i_species < n_species; i_species++){
        datafile << "," << species[i_species];
    }
    datafile << "\n";

    for (int i_cell = 0; i_cell < n_cells; i_cell++){
        datafile << i_cell;

        for (int i_species = 0; i_species < n_species; i_species++){
            datafile << "," << cells[i_cell].state[i_species];
        }

        datafile << "\n";
    }

    datafile.close();
}


void Model::error_dump(){

    cout << "Error dump generated" << endl;

    Reaction *reaction_ptr;

    string file_name;
    file_name = project_name + "_error_dump.txt";

    ofstream datafile;
    datafile.open(file_name);

    datafile << "All events and their propensities at time of error in " << model_name << " run " << i_run << endl;

    
    for(int i_event = 0; i_event < n_events; i_event++){
        
        reaction_ptr = events[i_event].reaction_ptr;
        datafile << i_event << " " << reaction_ptr->reaction_name << " in cell " << events[i_event].target_cells[0];
        if (reaction_ptr->is_diffusion_like){
            datafile << " to cell " << events[i_event].target_cells[1];
        }
        datafile << ", a = " << propensities[i_event] << "   " << endl;
    }

    datafile << "a0 = "<< a0 << endl;
    datafile << "Mu = "<< mu << endl;
    datafile << "Random number = "<< r_mu << endl;


	datafile << endl;
	datafile.close();

    // Also record concentrations normally
    record_concentrations(i_iter);
}

