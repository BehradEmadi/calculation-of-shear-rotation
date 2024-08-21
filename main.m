function shearrotation(nbeams, beamlength, measur_values, KnownName1, valuevars1)
    % Function to calculate shear rotation nbeams = number of beams, measur_values = deformation of each point knownname = beams with observed materail properties, valuevars= observed materail properties 
    global subbeamsEAIQI v % V = poisson value , subbeamsEAIQI = value of each element
    k = 5/6;  % Shear correction factor
    beamnode = nbeams + 1;
    nsize = 3 * beamnode;
    varnames = ['u'; 'v'; 'w'];
    allGSMVarNames = cell(nsize, 1);
    shearvars = cell(2 * nbeams, 1);
    count = 1;        
    cont2 = 0;

    % Generate shear variable names
    for j = 1:nbeams
        for beamn = 1:2
            cont2 = cont2 + 1;
            shearvars{cont2} = ['w' sprintf('%04.0f', count)];
            count = count + 1;
        end
        count = count - 1;
    end

    % Generate GSM variable names % GSM = global stiffness matric
    cont = 0;
    for k = 1:3
        for beamn = 1:beamnode
            cont = cont + 1;
            allGSMVarNames{cont} = [varnames(k,1) sprintf('%04.0f', beamn)];
        end
    end

    % Initialize measurements with zeros
    measurements = {allGSMVarNames, zeros(nsize, 1)};
    mearow = length(measur_values{2});

    % Assign measured values to the corresponding variables
    for ime = 1:mearow
        for jme = 1:nsize
            if strcmp(measur_values{1,1}(ime,1), measurements{1,1}(jme,1))
                measurements{1,2}(jme,1) = measur_values{1,2}(ime,1);
            end
        end
    end

    % Calculate shear force values
    [shearforcevalue, shearforcevars] = shearforce(beamnode, nbeams);

    % Initialize variables for known values
    auxE = []; auxI = []; auxA = [];
    valueE = []; valueI = []; valueA = [];
    cont1 = 1; cont2 = 1; cont3 = 1; cont4 = 1;

    % Separate known names into E, I, A categories
    for i = 1:length(KnownName1)
        aux1 = char(KnownName1(i));
        auxL = aux1(1:1);
        auxN = str2double(aux1(2:end));
        switch auxL
            case 'E'
                auxE(cont1) = auxN;
                valueE(cont1) = valuevars1(i);
                cont1 = cont1 + 1;
            case 'I'
                auxI(cont2) = auxN;
                valueI(cont2) = valuevars1(i);
                cont2 = cont2 + 1;
            case 'A'
                auxA(cont3) = auxN;
                valueA(cont3) = valuevars1(i);
                cont3 = cont3 + 1;
        end
    end

    % Calculate shear rotation and shear rotation by shear force
    SRBSF = zeros(2 * nbeams, 1);
    wr = zeros(nbeams, 1);

    for i = 1:length(subbeamsEAIQI{1})
        auxl = eval(beamlength{i});
        subE = subbeamsEAIQI{1}(i);
        subI = subbeamsEAIQI{3}(i);
        subA = subbeamsEAIQI{2}(i);
        v1 = measurements{1,2}(beamnode + i, 1);
        v2 = measurements{1,2}(beamnode + i + 1, 1);
        w1 = measurements{1,2}(2 * beamnode + i, 1);
        w2 = measurements{1,2}(2 * beamnode + i + 1, 1);

        % Find corresponding E, I, A values
        valE = valueE(auxE == subE);
        valI = valueI(auxI == subI);
        valA = valueA(auxA == subA);

        G = valE / (2 * (1 + v));
        a = (valE * valI) / auxl^3;
        b = k * G * valA * auxl;

        wr(i) = (-(12 * a * auxl) * v1 - (6 * a * auxl^2) * w1 + (12 * a * auxl) * v2 - (6 * a * auxl^2) * w2) ...
                / ((12 * a * auxl^2) + b);

        SRBSF(cont4) = shearforcevalue(cont4) / (k * valA * G); %SRBSF = shear rotation by shear force
        cont4 = cont4 + 1;
        SRBSF(cont4) = shearforcevalue(cont4) / (k * valA * G);
        cont4 = cont4 + 1;
    end

    % Evaluate wr and SRBSF wr= shear rotation
    wr = eval(wr);
    SRBSF = eval(SRBSF);

    % Transform shear rotation parameter to shear rotation
    shearresult = cell(2 * nbeams, 1);
    j = 1;
    for i = 1:nbeams
        shearresult{j} = wr(i);
        j = j + 1;
        shearresult{j} = wr(i);
        j = j + 1;
    end

    % Average shear rotation results
    avercont = zeros(beamnode, 1);
    shearvarsorg = unique(shearvars);
    shearresultorg = zeros(length(shearvarsorg), 1);
    for i = 1:length(shearvarsorg)
        for j = 1:length(shearvars)
            if isequal(shearvarsorg(i), shearvars(j))
                shearresultorg(i) = shearresultorg(i) + shearresult{j};
                avercont(i) = avercont(i) + 1;
            end
        end
    end
    shearresultorg = shearresultorg ./ avercont;

    % Transform shear rotation parameter by shear force to shear rotation
    SRV = unique(shearforcevars);
    SRRBSF = zeros(length(SRV), 1);
    avercont2 = zeros(beamnode, 1);

    for i = 1:length(SRRBSF)
        for j = 1:length(shearforcevars)
            if isequal(SRV(i), shearforcevars(j))
                SRRBSF(i) = SRRBSF(i) + SRBSF(j);
                avercont2(i) = avercont2(i) + 1;
            end
        end
    end
    SRRBSF = SRRBSF ./ avercont2;

    % Print results
    fprintf('\n--- Flexural Rotations\n');
    aaaa = measurements{1,1};
    for i = 1:beamnode
        fprintf('%s = %.10d\n', aaaa{2 * beamnode + i}, measurements{1,2}(2 * beamnode + i, 1));
    end

    fprintf('\n--- Shear Rotations\n');
    for i = 1:beamnode
        fprintf('%s = %.10d\n', aaaa{2 * beamnode + i}, shearresultorg(i));
    end

    fprintf('\n--- Shear Rotations by Shear force\n');
    for i = 1:beamnode
        fprintf('%s = %.10d\n', aaaa{2 * beamnode + i}, SRRBSF(i));
    end

    % Update measurements with shear results
    for i = 1:beamnode
        measurements{1,2}(2 * beamnode + i, 1) = measurements{1,2}(2 * beamnode + i, 1) + shearresultorg(i);
    end

    fprintf('\n--- Displacements\n');
    for i = 1:nsize
        fprintf('%s = %.10d\n', aaaa{i}, measurements{1,2}(i, 1));
    end
end
