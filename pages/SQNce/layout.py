import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, State, dash_table

PLOTLY_LOGO = "https://images.plot.ly/logo/new-branding/plotly-logomark.png"

layout = html.Div([ 
    dcc.Location(id='sqnce_url', refresh=False),
    # Use a row instead of a number and add cute buttons
    dbc.Navbar([
    #dbc.Container([
        html.Img(src=PLOTLY_LOGO, height="40px"),
        html.P("...", style={"color": "transparent"}),
        #dbc.NavItem(dbc.NavbarBrand("PlantApp")),
        #dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
        
        dbc.Button("About", id="about", outline=True, color="light", className="m-1", external_link=True, href="/about/"),
        dbc.DropdownMenu(
            [dbc.DropdownMenuItem("Available DBs", id="available_dbs", external_link=True, href="available_dbs"), 
            dbc.DropdownMenuItem(divider=True),
            dbc.DropdownMenuItem("Annotations", id="annotations", external_link=True, href="annotations"),
            dbc.DropdownMenuItem(divider=True),
            dbc.DropdownMenuItem("Gene Families", header=True),
            dbc.DropdownMenuItem("Using Family Names", id="families_familyIDs", external_link=True, href="families_familyIDs"),
            dbc.DropdownMenuItem("Using Gene IDs", id="families_geneIDs", external_link=True, href="families_geneIDs"),
            dbc.DropdownMenuItem(divider=True),
            dbc.DropdownMenuItem("Coordinates", id="coordinates", external_link=True, href="coordinates"),
            dbc.DropdownMenuItem(divider=True),
            dbc.DropdownMenuItem("Gene Sequences", header=True),
            dbc.DropdownMenuItem("Sequences", id="proteins", external_link=True, href="proteins"),
            dbc.DropdownMenuItem("Promoters", id="promoters", external_link=True, href="promoters"),
            dbc.DropdownMenuItem(divider=True),
            dbc.DropdownMenuItem("Transcriptomic Data", header=True),
            dbc.DropdownMenuItem("Omics", id="omics", external_link=True, href="omics"),
            ], label="SQNce", color="primary", className="m-1",
        ),
        dbc.DropdownMenu(
            [dbc.DropdownMenuItem("Simple Tree", id="simple_tree", external_link=True, href="/apps/"),
            dbc.DropdownMenuItem("GO Enrichment", id="go_enrichment", disabled=True),
            ], label="Apps", color="secondary", className="m-1",
            toggle_style={"background": "transparent","border-color": "#f8f9fa"},
        ),
        dbc.Button("Downloads", id="downloads", outline=True, color="secondary", 
                    className="m-1", disabled=True),
        dbc.Button("Feedback", id="feedback", outline=True, color="secondary", 
                    className="m-1", disabled=True),
        ],
        #align="start",
        #style={"margin-left": "15px"}
        dark =True, 
        color="dark",
        sticky = True,
    ),
    html.Div(id='page-content'),
    html.Div(id="tab_content"),
    ],
    #fluid=True,
)