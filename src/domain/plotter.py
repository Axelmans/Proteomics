import matplotlib.pyplot as plot
import sqlite3
from src.domain.values import EthnicGroups

CONN = sqlite3.connect("src/db/proteomics.db")

# This module generates graphs for the PTM match frequencies of individuals.
# It is expected that you run this script after filling the database with PTM matches!
if __name__ == '__main__':
    CURSOR = CONN.cursor()
    # Iterate over all individual numbers and plot a graph for each.
    for group in EthnicGroups:
        for number in group.value:
            print(f"Analysing individual {number}.")
            # The match frequencies are queried from the database.
            CURSOR.execute(
                "SELECT * FROM ptm_match_frequencies WHERE individual = ?",
                (number,)
            )
            items = CURSOR.fetchall()
            # Each match and its frequency gets an entry in the plot.
            ptms = []
            frequencies = []
            for item in items:
                ptms.append(item[2])
                frequencies.append(item[3])
            plot.figure(figsize=(10, 6))
            plot.bar(ptms, frequencies)
            # Ptms go on the x-axis, frequencies on the y-axis.
            plot.xlabel("PTM")
            plot.ylabel("Frequency")
            plot.title(f"PTM Match Frequencies for Individual {number} ({group.name.lower()})")
            plot.xticks(rotation=45, ha='right')
            plot.tight_layout()
            # Plots are saved as .png files in the ptm_match_graphs directory
            plot.savefig(f"ptm_match_graphs/{group.name.lower()}/{number}.png")
            plot.close()
            print("Done!")
