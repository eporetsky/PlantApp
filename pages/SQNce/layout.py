import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, State, dash_table

PLOTLY_LOGO = "https://images.plot.ly/logo/new-branding/plotly-logomark.png"

layout = dbc.Container(
    [
        dcc.Store(id="store"),

        #html.H1("Dynamically rendered tab content"),
        #html.Hr(),

    dbc.Navbar(
        children = [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [dbc.Col(html.Img(src=PLOTLY_LOGO, height="30px")),
                    dbc.Col(dbc.NavbarBrand("PlantApp", className="ml-2")),
                    dbc.Col(dbc.NavLink("About", href="/about/", className="ml-2", external_link=True)),
                    dbc.Col(dbc.NavLink("SQNce", href="/SQNce/", className="ml-2", external_link=True)),
                    dbc.Col(dbc.NavLink("Queries", href="/queries/", className="ml-2", external_link=True)),
                    dbc.Col(dbc.NavLink("Tools", href="/tools/", className="ml-2", external_link=True)),
                    dbc.Col(dbc.NavLink("Downloads", href="/downloads/", className="ml-2", external_link=True)),
                    dbc.Col(dbc.NavLink("Feedback", href="/feedback/", className="ml-2", external_link=True)),
                    ],
                    align="center",# no_gutters=True,
                ),
                href="https://www.plantapp.org",
            ),
            dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
            #dbc.Collapse(search_bar, id="navbar-collapse", navbar=True, is_open=False),
        ],
    color="dark", dark=True,
    ),
    
    dbc.Tabs(
        [
            dbc.Tab(label="Available DBs", tab_id="available_dbs"),
            dbc.Tab(label="Annotations", tab_id="annotations"),
            dbc.Tab(label="Coordinates", tab_id="coordinates"),
            dbc.Tab(label="Sequences", tab_id="proteins"),
            dbc.Tab(label="Promoters", tab_id="promoters"),
            dbc.Tab(label="Omics", tab_id="omics"),
        ],
        id="tabs",
        active_tab="scatter",
    ),
        
    html.Div(id="tab_content"),

        
      #dbc.Row([table])
    ],
    fluid=True  
)