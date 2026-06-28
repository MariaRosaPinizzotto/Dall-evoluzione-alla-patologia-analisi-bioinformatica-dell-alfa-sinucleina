%% Progetto Bioinformatica: Analisi dell'Alfa-Sinucleina
% Passo 1: Caricamento dati (Versione Compatibile)

clear;       % Cancella le vecchie variabili
clc;         % Pulisce lo schermo sotto

fprintf('Inizio analisi Alfa-Sinucleina...\n');

% Percorso del file
file_path = fullfile('..', 'data', 'alfa_sinucleina_human.fasta');

% Leggiamo il file riga per riga in modo standard
fid = fopen(file_path, 'r');
if fid == -1
    error('Accidenti! Il file FASTA non è stato trovato nella cartella data.');
end

header = fgetl(fid); % Legge la prima riga (l'intestazione)
sequenza = '';

% Legge tutte le altre righe attaccandole insieme
while ~feof(fid)
    linea = strtrim(fgetl(fid));
    sequenza = [sequenza, linea];
end
fclose(fid);

% Stampare a video i dati letti
fprintf('Proteina caricata con successo!\n');
fprintf('Intestazione: %s\n', header);
fprintf('Numero di amminoacidi: %d\n', length(sequenza));

%% Passo 2: Analisi della composizione amminoacidica
fprintf('\n--- Analisi della Composizione ---\n');

% Trovare tutti gli amminoacidi unici presenti e ordinarli
amino_acids = unique(sequenza);
conteggi = zeros(length(amino_acids), 1);

% Numerazione di amminoaicidi per tipo
for i = 1:length(amino_acids)
    conteggi(i) = sum(sequenza == amino_acids(i));
end

% Stampare a video i risultati nella Command Window
disp('Amminoacido -> Quantità:');
for i = 1:length(amino_acids)
    fprintf('   %s -> %d\n', amino_acids(i), conteggi(i));
end

% Creazione il grafico a barre per visualizzare il risultato
figure('Name', 'Composizione Amminoacidica Alfa-Sinucleina');
bar(conteggi, 'FaceColor', [0 0.4470 0.7410]);
set(gca, 'XTick', 1:length(amino_acids), 'XTickLabel', cellstr(amino_acids'));
title('Composizione Amminoacidica dell''Alfa-Sinucleina Umana');
xlabel('Amminoacidi (Codice a una lettera)');
ylabel('Numero di occorrenze');
grid on;

% Salvataggio automatico del grafico nella cartella 'figures'
saveas(gcf, fullfile('..', 'figures', 'composizione_amminoacidica.png'));
fprintf('\nGrafico salvato correttamente nella cartella "figures"!\n');

%% NOTE BIOLOGICHE SUI DIVERSI AMMINOACIDI DELL'ALFA-SINUCLEINA

% 1. AMMINOACIDI CARICHI POSITIVAMENTE - LISINA (K)
% La Lisina (K) è molto frequente nella regione N-terminale. 
% Le sue cariche positive permettono alla proteina di legarsi 
% ai fosfolipidi (carichi negativamente) delle membrane delle vescicole sinaptiche.

% 2. AMMINOACIDI PICCOLI E IDROFOBICI - ALANINA (A), GLICINA (G), VALINA (V)
% Alanina e Glicina sono abbondantissime. La loro piccola dimensione dona 
% grande flessibilità alla proteina. Nella regione centrale (chiamata NAC), 
% amminoacidi idrofobici come la Valina (V) tendono a "scappare" dall'acqua; 
% quando la proteina si ripiega male, è proprio questa zona che si incolla 
% ad altre molecole creando gli aggregati patologici (corpi di Lewy).

% 3. AMMINOACIDI CARICHI NEGATIVAMENTE - ACIDO GLUTAMMICO (E)
% L'Acido Glutammico (E) è concentrato nella regione C-terminale (la coda). 
% Essendo carico negativamente, crea repulsione elettrostatica che impedisce 
% alla proteina di aggrovigliarsi su se stessa in condizioni normali.

% 4. AMMINOACIDI ASSENTI O RARI - TRIPTOFANO (W) E CISTEINA (C)
% Nota che nel grafico la Cisteina (C) e il Triptofano (W) sono assenti (0). 
% L'assenza di Cisteine significa che l'alfa-sinucleina non può formare 
% ponti disolfuro stabili, il che spiega perché rimane una proteina disordinata.