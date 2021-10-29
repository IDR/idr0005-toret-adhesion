#! /usr/bin/env python

import csv
import logging
import os
DEBUG = int(os.environ.get("DEBUG", logging.INFO))
logging.basicConfig(level=DEBUG)


def read_gene_mappings(screen):
    """Parse a validation table from Flybase into a dictionary"""
    genes = {}
    nmappings = 0
    with open('%s/gene_mappings.txt' % screen, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        next(reader, None)  # skip the headers
        for row in reader:
            nmappings += 1
            # There could be multiple target genes associated to an entry
            genes.setdefault(row[0], set()).add((row[1], row[2]))
    logging.info("Found %s validated genes (%s hits in total)" %
                 (len(genes), nmappings))
    return genes


def find_match(candidates, symbol):
    # If several candidate exist, return the first one matching the gene symbol
    if (len(candidates) > 1):
        for candidate in candidates:
            if (candidate[1] == symbol):
                return candidate
        return
    else:
        return next(iter(candidates))


def normalise_genes(screen):
    n_genes = 0
    n_updates = 0
    n_symbol_updates = 0
    n_synonyms_updates = 0
    genes = read_gene_mappings(screen)
    with open('%s/idr0005-%s-annotation.csv' % (screen, screen),
              newline='', mode='r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        headers = next(csv_reader, None)
        # Add new columns for validation
        headers.insert(16, 'New Gene Symbol')
        headers.insert(17, 'New Gene Synonyms')
        data = []
        for row in csv_reader:
            row.insert(16, '')
            row.insert(17, '')
            gene_id = row[13]
            if gene_id.startswith('FBgn') or gene_id not in genes:
                data.append(row)
                continue

            n_genes += 1
            gene_symbol = row[14]
            gene_synonyms = row[15]

            match = find_match(genes[gene_id], gene_symbol)
            if not match:
                data.append(row)
                continue

            new_gene_id = match[0]
            new_gene_symbol = match[1]

            if new_gene_symbol == 'unknown ID':
                data.append(row)
                continue

            # If the validated symbol has changed
            if (new_gene_symbol != gene_symbol):
                row[16] = new_gene_symbol
                n_symbol_updates += 1

            if gene_id not in gene_synonyms:
                # If the old identifier if not in the list of synonyms
                if gene_synonyms == '':
                    row[17] = gene_id
                else:
                    row[17] = gene_synonyms + ' ' + gene_id
                n_synonyms_updates += 1
            row[13] = new_gene_id
            data.append(row)
            n_updates += 1

    with open('%s/idr0005-%s-annotation-validated.csv' % (screen, screen),
              newline='', mode='w') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        csv_writer.writerow(headers)
        csv_writer.writerows(data)
        logging.info("Updated %s/%s gene rows with %s gene symbol updates "
                     "and %s gene synonym updates" %
                     (n_updates, n_genes, n_symbol_updates,
                      n_synonyms_updates))


normalise_genes('screenA')
normalise_genes('screenB')
