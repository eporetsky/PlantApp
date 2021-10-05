import dash_bootstrap_components as dbc

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