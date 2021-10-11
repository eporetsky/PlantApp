import dash_bootstrap_components as dbc

def dropdown_menus():
    con = sqlite3.connect('SQNce.db')
    cursorObj = con.cursor()
    # read in your SQL query results using pandas
    dropdown_df = pd.read_sql("""
        SELECT  studies.study_accession, studies.tax_id, fastq.run_accession, fastq.sample_alias
        FROM studies, fastq
        WHERE studies.study_accession=fastq.study_accession and studies.tax_id==4577
        """, con = con)
    return(dropdown_df)



search_bar = [ 
    dbc.Row([
    dbc.Col(dbc.Input(type="search", placeholder="Search")),
        dbc.Col(
            dbc.Button("Search", color="primary", className="ml-2", n_clicks=0),
            width="auto",),
        ],
    no_gutters=True,
    className="ml-auto flex-nowrap mt-3 mt-md-0",
    align="center",
    ), 
]

