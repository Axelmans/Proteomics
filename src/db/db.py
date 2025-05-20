import sqlite3

CONN = sqlite3.connect("src/db/proteomics.db")
CURSOR = CONN.cursor()

# Table creation.
CURSOR.execute(
    """
    CREATE TABLE IF NOT EXISTS ptm_match_frequencies (
        individual INTEGER,
        ethnic_group TEXT,
        ptm_match TEXT,
        frequency INTEGER,
        PRIMARY KEY (individual, ptm_match)
    )
    """
)

def clear_tables():
    CURSOR.execute("DELETE FROM ptm_match_frequencies")

# Adding a match does not require an input frequency, this value is updated automatically.
def add_ptm_match(individual, ethnic_group, ptm_match):
    CURSOR.execute(
    """
        INSERT INTO ptm_match_frequencies (individual, ethnic_group, ptm_match, frequency)
        VALUES (?, ?, ?, 1)
        ON CONFLICT(individual, ptm_match)
        DO UPDATE SET frequency = frequency + 1
        """,(individual, ethnic_group, ptm_match)
    )
    CONN.commit()

# Test adding entries
if __name__ == "__main__":
    clear_tables()
    add_ptm_match(1, "White", "METHYLATION")
    add_ptm_match(1, "White", "METHYLATION")
