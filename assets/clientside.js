// assets/clientside.js
window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        toggleTheme: function(data) {
            const theme = data.theme || 'light';
            document.body.className = theme + "-theme";
            console.log('Theme switched to', theme);
            return '';
        }
    }
});
