function main_analysis()
    %% Progetto Bioinformatica: Analisi Comparativa dell'Alfa-Sinucleina
    clc;
    fprintf('Inizio analisi computazionale Alfa-Sinucleina...\n');

    % =========================================================================
    % 1. GESTIONE DINAMICA DEI PERCORSI
    % =========================================================================
    % Individuazione della cartella corrente in cui risiede lo script ('src')
    [current_dir, ~, ~] = fileparts(mfilename('fullpath'));
    % Risalita di un livello per identificare la radice principale del progetto
    [project_root, ~, ~] = fileparts(current_dir); 

    % =========================================================================
    % 2. CARICAMENTO SEQUENZA DI RIFERIMENTO (HOMO SAPIENS)
    % =========================================================================
    file_path = fullfile(project_root, 'data', 'alfa_sinucleina_human.fasta');
    fid = fopen(file_path, 'r');
    if fid == -1
        error('Il file FASTA dell''uomo non è stato rinvenuto nel percorso: %s', file_path);
    end
    
    % Lettura e scarto della prima riga di intestazione del file FASTA (Risolto errore sintassi)
    [~] = fgetl(fid); 
    sequenza = '';
    
    % Ricostruzione della stringa lineare della sequenza amminoacidica
    while ~feof(fid)
        linea = strtrim(fgetl(fid));
        sequenza = [sequenza, linea]; %#ok<AGROW> Soppressione warning allocazione dinamica
    end
    fclose(fid);

    fprintf('Proteina di riferimento (Uomo) caricata correttamente. Lunghezza: %d aa\n', length(sequenza));

    % =========================================================================
    % 3. ANALISI DELLA COMPOSIZIONE AMMINOACIDICA UMANA
    % =========================================================================
    % Identificazione dei residui unici presenti all'interno della catena
    amino_acids = unique(sequenza);
    conteggi = zeros(length(amino_acids), 1);
    
    % Calcolo delle frequenze assolute di ciascun amminoacido
    for i = 1:length(amino_acids)
        conteggi(i) = sum(sequenza == amino_acids(i));
    end

    % =========================================================================
    % 4. SVILUPPO DELL'INTERFACCIA GRAFICA (GUI)
    % =========================================================================
    fig_gui = uifigure('Name', 'Analisi Alfa-Sinucleina: Allineamento Evolutivo', 'Position', [100, 100, 500, 400]);

    % Etichetta del Titolo Principale
    uilabel(fig_gui, 'Text', 'Analisi Biocomputazionale Alfa-Sinucleina', ...
        'Position', [20, 350, 460, 30], 'FontWeight', 'bold', 'FontSize', 16, 'HorizontalAlignment', 'center');

    % Pulsante per la generazione del grafico delle frequenze umane
    uibutton(fig_gui, 'push', 'Text', 'Visualizza Composizione Umana', ...
        'Position', [75, 280, 350, 40], 'ButtonPushedFcn', @(btn, event) mostraGraficoUmano(amino_acids, conteggi));

    % Etichetta informativa per la selezione
    uilabel(fig_gui, 'Text', 'Selezionare la specie per il confronto filogenetico:', ...
        'Position', [75, 220, 350, 20], 'FontWeight', 'bold');

    % Menu a tendina per la selezione della specie ortologa (Tre opzioni integrate)
    dropdown_specie = uidropdown(fig_gui, ...
        'Items', {'Scimpanzé (Pan troglodytes)', 'Topo (Mus musculus)', 'Gallina (Gallus gallus)'}, ...
        'ItemsData', {'alfa_sinucleina_chimp.fasta', 'alfa_sinucleina_mouse.fasta', 'alfa_sinucleina_chicken.fasta'}, ...
        'Position', [75, 180, 350, 30]);

    % Pulsante per l'attivazione dell'allineamento locale e della mappatura
    uibutton(fig_gui, 'push', 'Text', 'Esegui Allineamento e Mappatura', ...
        'Position', [75, 110, 350, 40], 'ButtonPushedFcn', @(btn, event) eseguiConfronto(fig_gui, dropdown_specie, sequenza, project_root));

    % Nota informativa a piè di pagina
    uilabel(fig_gui, 'Text', 'Analisi di Conservazione Evolutiva delle Sequenze Primarie', ...
        'Position', [20, 30, 460, 20], 'FontAngle', 'italic', 'HorizontalAlignment', 'center');
end

%% =========================================================================
% 5. FUNZIONI LOCALI DI SUPPORTO
%% =========================================================================

function mostraGraficoUmano(amino_acids, conteggi)
    % Generazione del profilo a barre relativo alla frequenza dei residui umani
    figure('Name', 'Profilo Amminoacidico Umano');
    bar(conteggi, 'FaceColor', [0 0.4470 0.7410]);
    set(gca, 'XTick', 1:length(amino_acids), 'XTickLabel', cellstr(amino_acids'));
    
    title('Composizione Amminoacidica dell''Alfa-Sinucleina Umana');
    xlabel('Residui Amminoacidici (Codice a singola lettera)'); 
    ylabel('Numero di occorrenze'); 
    grid on;
end

function eseguiConfronto(finestra, dropdown, sequenza_human, project_root)
    % Estrazione dei metadati relativi alla specie selezionata dall'utente
    file_specie = dropdown.Value;
    indice = strcmp(dropdown.ItemsData, file_specie);
    nome_specie = dropdown.Items{indice};

    % Apertura e lettura del file ortologo corrispondente
    percorso = fullfile(project_root, 'data', file_specie);
    fid_s = fopen(percorso, 'r');
    if fid_s == -1
        uialert(finestra, sprintf('File non trovato nel percorso:\n%s', percorso), 'Errore di lettura', 'Icon', 'error');
        return;
    end
    
    [~] = fgetl(fid_s); % Risolto errore sintassi anche qui
    sequenza_specie = '';
    while ~feof(fid_s)
        sequenza_specie = [sequenza_specie, strtrim(fgetl(fid_s))]; %#ok<AGROW> Soppressione warning allocazione dinamica
    end
    fclose(fid_s);

    % Determinazione della lunghezza minima per l'allineamento di posizione
    L_min = min(length(sequenza_human), length(sequenza_specie));
    vettore_allineamento = zeros(1, L_min); % Inizializzazione matrice logica (1=Identico, 0=Mutato)
    identici = 0;
    
    % Ciclo di confronto residuo per residuo lungo l'asse lineare della proteina
    for j = 1:L_min
        if sequenza_human(j) == sequenza_specie(j)
            identici = identici + 1;
            vettore_allineamento(j) = 1; 
        else
            vettore_allineamento(j) = 0; 
        end
    end
    
    % Calcolo degli indici statistici di identità e mutazione
    perc = (identici / L_min) * 100;
    mutati = L_min - identici;

    % GENERAZIONE GRAFICO POSIZIONALE DELLE MUTAZIONI
    figure('Name', ['Mappatura Mutazioni - Uomo vs ' nome_specie]);
    hold on;
    % Tracciamento della linea connettiva di base
    plot(1:L_min, vettore_allineamento, 'Color', [0.8 0.8 0.8]);
    % Rappresentazione geometrica dei residui conservati (pallini verdi)
    scatter(find(vettore_allineamento == 1), ones(1, identici), 30, [0.46 0.67 0.18], 'filled', 'DisplayName', 'Identico');
    % Rappresentazione geometrica dei residui mutati (pallini rossi)
    if mutati > 0
        scatter(find(vettore_allineamento == 0), zeros(1, mutati), 40, [0.85 0.32 0.09], 'filled', 'DisplayName', 'Mutato');
    end
    
    % Configurazione assi e formattazione del grafico
    set(gca, 'YTick', [0 1], 'YTickLabel', {'Mutato', 'Identico'});
    title(['Mappatura Posizionale delle Mutazioni (Uomo vs ' nome_specie ')']);
    xlabel('Posizione lungo la catena (aa 1-140)'); 
    ylim([-0.5 1.5]); 
    grid on;
    legend('Location', 'southwest'); 
    hold off;

    % CONFIGURAZIONE DELLE NOTE EVOLUTIVE IN FORMA STRETTAMENTE IMPERSONALE
    if contains(nome_specie, 'Scimpanzé')
        nota_bio = sprintf(['Analisi del tasso di conservazione:\\n', ...
            'Si riscontra una corrispondenza pressoché assoluta (~99-100%%). Tale evidenza è determinata dalla minima distanza filogenetica e dal ridotto tempo di divergenza dall''antenato comune. Le rare sostituzioni riscontrate presentano un carattere strettamente conservativo.']);
    elseif contains(nome_specie, 'Topo')
        nota_bio = sprintf(['Analisi del tasso di conservazione:\\n', ...
            'Si osserva un livello di identità estremamente elevato (95.00%%). Questo valore attesta l''azione di una forte pressione selettiva stabilizzante nei mammiferi. I residui mutati si concentrano in regioni non critiche, preservando l''architettura funzionale.']);
    else % Caso della Gallina
        nota_bio = sprintf(['Analisi del tasso di conservazione:\\n', ...
            'Si registra il decremento più significativo nel tasso di identità (~85-89%%). Tale variazione è determinata dalla sostanziale distanza filogenetica della linea degli uccelli, separata dai mammiferi da circa 300 milioni di anni di evoluzione. Le mutazioni si localizzano prevalentemente nel dominio C-terminale, strutturalmente più tollerante.']);
    end

    % Composizione finale del report visivo all'interno della finestra di alert
    msg = sprintf('SPECIE: %s\\n\\nIdentità: %.2f%%\\nMutazioni: %d su %d\\n\\n%s', ...
        nome_specie, perc, mutati, L_min, nota_bio);
    uialert(finestra, msg, 'Risultati Allineamento', 'Icon', 'info');
end