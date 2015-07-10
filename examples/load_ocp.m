function p = load_ocp(filename)

p = struct('n_states', 0, 'n_controls', 0, ...
    'n_parameters', 0, 'n_initial', 0, 'n_terminal', 0, ...
    'n_inequality', 0, 'tolerance', 0, 'c_type', 2,
    'Y', [], 'U', [], 'P', [], ...
    'T', []);

    file_id = fopen(filename, 'r');
    
    s = fscanf(file_id, '%s',1);
    i = fscanf(file_id, '%d',1);
    p.n_states = i;
    s = fscanf(file_id, '%s',1);
    i = fscanf(file_id, '%d',1);
    p.n_controls = i;
    s = fscanf(file_id, '%s',1);
    i = fscanf(file_id, '%d',1);
    p.n_parameters = i;
    s = fscanf(file_id, '%s',1);
    i = fscanf(file_id, '%d',1);
    p.n_initial = i;
    s = fscanf(file_id, '%s',1);
    i = fscanf(file_id, '%d',1);
    p.n_terminal = i;
    s = fscanf(file_id, '%s',1);
    i = fscanf(file_id, '%d',1);
    p.n_inequality = i;
    s = fscanf(file_id, '%s',1);
    i = fscanf(file_id, '%d',1);
    p.n_nodes = i;
    s = fscanf(file_id, '%s',1);
    x = fscanf(file_id, '%le',1);
    p.tolerance = x;
    s = fscanf(file_id, '%s',1);
    i = fscanf(file_id, '%d',1);
    p.c_type = i;
    
    p.T = zeros(p.n_nodes, 1);
    p.Y = zeros(p.n_states, p.n_nodes);
    if(p.n_parameters>0)
        p.P = zeros(p.n_parameters, 1);
    end;
    if(p.n_controls>0)
        p.U = zeros(p.n_controls, p.n_nodes);
    end;
    s = fscanf(file_id, '%s',1); % /* T */
    for i=1:p.n_nodes
        x = fscanf(file_id, '%le',1);
        p.T(i) = x;
    end;
    s = fscanf(file_id, '%s',1); % /* Y */
    for i=1:p.n_states
        for j=1:p.n_nodes
            x = fscanf(file_id, '%le',1);
            p.Y(i,j) = x;
        end;
    end;
    s = fscanf(file_id, '%s',1); %/* U */
    for i=1:p.n_controls
        for j=1:p.n_nodes
            x = fscanf(file_id, '%le',1);
            p.U(i,j) = x;
        end;
    end;
    s = fscanf(file_id, '%s',1); %/* P */
    for i=1:p.n_parameters
        x = fscanf(file_id, '%le',1);
        p.P(i) = x;
    end;


